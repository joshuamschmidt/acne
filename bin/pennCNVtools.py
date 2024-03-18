#!/usr/bin/env python
import argparse
import numpy as np
import polars as pl
import textwrap
import os
import re

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='pennCNVtools',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent('''\
                        Python Utilities for making PENN CNV inputs from
                        Genome Studio output (GType, BAF and Log R Ratio)
                        ------------------------------------------------
                        Option "partition" splits into partitions. Useful
                        for large input files (n ind > 500). It will try
                        to balance e.g. if n ind = 560, and n_per_partition
                        is 250 it will produce two files with 280 inds,
                        not 250 + 310 nor 250 + 250 + 60. 


                        Option "pfb" Generates a PFB file for PENN CNV

                        Output Column headers are:
                        Name    Chr     Position        PFB 
                        
                        Option "split" Generates a ind file for
                        PENN CNV

                        Output Column headers are:
                        Name  Sample.GT sample.Log R Ratio Sample.BAF

                        Note - these tools use Polars for fast data input
                        and manipulation. set env for n processes.
                        '''),
                        add_help=False,
                        epilog="Questions, bugs etc?\njoshmschmidt1@gmail.com\ngithub.com/joshuamschmidt")
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

# Add back help
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=argparse.SUPPRESS,
    help='show this help message and exit'
)

parser.add_argument('tool',metavar='TOOL', type=str, nargs=1, choices={"partition", "split", "pfb"}, help='which tool to run')

required.add_argument('--input', type=str, dest='input',
                    help='input file from genome studio with GT, LogR and BAF cols per sample')

optional.add_argument('--output', type=str,
                       dest='output',
                       help='output:\nif tool is "pfb", then creates a pfb file\nit has no effect for other tools')


optional.add_argument('--n', type=int, dest='n_per_partition',
                      help='how many inds should large input be split into',
                      default=None)

optional.add_argument('--prefix', type=str, dest='prefix',
                      help='batch prefix to attach to ind sample files when splitting',
                      default=None)

optional.add_argument('--geno', type=float, dest='geno',
                      help='SNPs/markers with missingness fraction more than this are excluded from PFB',
                      default=0.02)

#---- functions used by multiple classes....

def file_str(input):
    file_str = {}
    with open(input, 'rt') as fh:
        header=fh.readline().strip()
    assert 'Name' in header, 'BAF file must have SNP Name column'
    std_cols = ['Name', 'Chr', 'Pos']
    file_str['n_std']=0
    file_str['std_cols']=[]
    for c in std_cols:
         match = re.search(c+'\t|'+c+',', header)
         if match:
            file_str['n_std'] += 1
            file_str['std_cols'].append(c)
    file_str['n_BAF']=header.count('B Allele Freq')
    file_str['n_LRR']=header.count('Log R Ratio')
    file_str['n_GT']=header.count('GType')
    assert file_str['n_BAF'] >= 1, 'You must have BAF data for at least 1 sample'
    assert file_str['n_LRR'] >= 1, 'You must have LRR data for at least 1 sample'
    assert file_str['n_BAF'] == file_str['n_LRR'], "n LRR and n BAF mismatch"
    assert file_str['n_GT'] == 0 or file_str['n_GT'] == file_str['n_BAF'], 'n GType and n BAF mismatch'
    file_str['exp_n'] = file_str['n_std'] + file_str['n_BAF'] + file_str['n_LRR'] + file_str['n_GT']
    return(file_str)
    



    


def make_clean_cols(input):
    
    cols=pl.read_csv(input, sep=sep_char, has_header=True, n_rows=1, n_threads=1, infer_schema_length=1 )
    if cols.shape[1]==1:
        sep_char='\t'
        cols=pl.read_csv(input, sep=sep_char, has_header=True, n_rows=1, n_threads=1, infer_schema_length=1 )
    cols=cols.columns
    clean_cols=[col if (col not in cols[:i]) else "DUP"+str(cols[:i].count(col))+"_"+str(col) for i, col in enumerate(cols)]
    return(clean_cols)

def get_col_types(input):
    col_types = {}




#----- class defs

# '''class for data to split into n ind chunks'''
class sampleDataPartition():
    def __init__(self, input: str, target_n: int):
        self.input = input
        self.target_n = target_n
        self.df = []
        self.clean_cols = []
        self._make_clean_cols()
        self.load_data()
        self.samples= []
        self.get_samples()
        self.n_samples = len(self.samples) # total number in input file
        self.n_per_partition=[]
        self.define_partition_n()
        self.prefix = os.path.splitext(input)[0]
    
    def _make_clean_cols(self.input):
        make_clean_cols(self)

    def load_data(self):
        q = (
            pl.scan_csv(self.input, separator='\t')
            .sort([
                pl.col("Chr"), pl.col("Position")],
                )
        )
        self.df = q.collect()
    
    def get_samples(self):
        self.samples=self.df.select(pl.col("^*.GType$")).columns
        self.samples=[s.split('.GType')[0] for s in self.samples]

    def define_partition_n(self):
        if(self.n_samples < 2 * self.target_n):
            print("target_n is larger than input sample_n / 2. Defaulting to splitting input in half")
            self.target_n, remainder = divmod(self.n_samples, 2)
            if(remainder > 0):
                self.target_n += 1
        n_partitions, remainder = divmod(self.n_samples, self.target_n)
        self.n_per_partition = [self.target_n] * n_partitions
        if(remainder > n_partitions):
            add_all, remainder = divmod(remainder, n_partitions)
            self.n_per_partition = [n + add_all for n in self.n_per_partition]
        if(remainder > 0):
            for i in range(remainder):
                self.n_per_partition[i] += 1

    def write_partition_data(self):
        all_samples = self.samples
        for j, n_part in enumerate(self.n_per_partition):
            part_samples = all_samples[:n_part]
            part_cols = [[s+".GType", s+".Log R Ratio", s+".B Allele Freq"] for s in part_samples]
            part_cols = [item for sublist in part_cols for item in sublist]
            sub=self.df.select(["Name", "Chr", "Position"] + part_cols)
            sub.write_csv(self.prefix + "-" + str(j+1) + '.partition', separator='\t')
            all_samples = all_samples[n_part:]




# '''class for data to split by ind'''
class sampleDataSplit():
    def __init__(self, input: str, prefix: str):
        self.input = input
        self.prefix = prefix
        self.clean_cols = []
        self._make_clean_cols()
        self.df = []
        self.load_data()
        self.samples= []
        self.get_samples()
    
    def _make_clean_cols(self):
        make_clean_cols(self.input)

    def load_data(self):
        q = (
            pl.scan_csv(self.input, separator='\t')
            .sort([
                pl.col("Chr"), pl.col("Position")],
                )
        )
        self.df = q.collect()
    
    def get_samples(self):
        self.samples=self.df.select(pl.col("^*.GType$")).columns
        self.samples=[s.split('.GType')[0] for s in self.samples]
    
    def write_sample_data(self):
        for s in self.samples:
            col_1, col_2, col_3 = s+".GType", s+".Log R Ratio", s+".B Allele Freq"
            sub=self.df.select(["Name","Chr","Position", col_1, col_2, col_3])
            sub.write_csv(self.prefix+'_'+s+'.txt', separator='\t')


# '''class for GtLogRBaf to pfb'''
class pfbObj():
    def __init__(self, input: str, geno: float,):
        self.input = input
        self.geno = geno
        self.clean_cols = []
        self._make_clean_cols()
        self.get_pfb()

    def _make_clean_cols(self):
        make_clean_cols(self.input)

    def get_pfb(self):
        col_types = [ pl.Utf8 if col=='Name' else pl.Categorical if col=='Chr' else pl.UInt32 if col=='Position' else pl.Categorical if 'GType' in col else pl.Float32 if 'Freq' in col else pl.Float32 if 'Ratio' in col else pl.UInt16 for col in self.clean_cols]
        col_dict = dict(zip(self.clean_cols, col_types))
        q = (
            pl.scan_csv(self.input, separator='\t', has_header=False, skip_rows=1, with_column_names=lambda cols: self.clean_cols, dtypes = col_dict)
            .select([pl.col("Name"), pl.col("Chr"), pl.col("Position"), pl.col("^*.B Allele Freq$")])
            .with_columns([
                pl.sum_horizontal(pl.col("^*.B Allele Freq$").is_nan()).alias('n_miss'),
                pl.sum_horizontal(pl.col("^*.B Allele Freq$").is_not_nan()).alias('n_call'),
                ])
            .fill_nan(0)
            .with_columns([
                pl.fold(acc=pl.lit(0),
                function=lambda acc, x: acc + x, exprs=pl.col("^*.B Allele Freq$")).alias("sum")
                ])
            .with_columns([
                ( (pl.col("sum") / (pl.col("n_call") - pl.col("n_miss")) )).alias("mean")
                ])
            .fill_nan(0)
            .with_columns([
                ((pl.col("mean")*1000+0.5).cast(pl.Int64)/1000).alias("BAF")
                ])
            .select(["Name","Chr","Position","BAF","n_miss","n_call"])
            .sort([
                pl.col("Chr"), pl.col("Position")],)
            )
        s= q.collect( streaming=True  )
        s=s.filter(pl.col("n_miss")/(pl.col("n_miss")+pl.col("n_call")) < self.geno)
        s=s.select([
            "Name",
            "Chr",
            "Position",
            pl.when(pl.col("Name").str.contains("cnv|CNV")).then(pl.lit(2)).otherwise(pl.col("BAF")).alias("PFB")
            ])
        self.pfb = s

def main():
    args = parser.parse_args()
    tool=args.tool[0]
    if(tool=='partition'):
        data = sampleDataPartition(args.input,args.n_per_partition)
        data.write_partition_data()

    if(tool=='split'):
        data = sampleDataSplit(args.input, args.prefix)
        data.write_sample_data()

    if(tool=='pfb'):
        if not args.output:
            parser.error('pfb tool selected: --output must be specified')
        pfb = pfbObj(args.input, args.geno)
        pfb.pfb.write_csv(args.output, separator='\t')


if __name__ == '__main__':
    main()
