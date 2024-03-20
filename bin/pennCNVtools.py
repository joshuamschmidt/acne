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

optional.add_argument('--map_file', type=str, dest='map_file',
                      help='file of Name, Chr, Pos to update co-ordinates in the PFB file',
                      default=None)


#---- functions used by multiple classes....

def make_file_struct(obj):
    file_struct = {}
    with open(obj.input, 'rt') as fh:
        file_struct['header']=fh.readline().strip()
    assert 'Name' in file_struct['header'], 'BAF file must have SNP Name column'
    std_cols = ['Name', 'Chr', 'Pos']
    file_struct['n_std']=0
    file_struct['std_cols']=[]
    for c in std_cols:
         match = re.search(c+'\t|'+c+',', file_struct['header'])
         if match:
            file_struct['n_std'] += 1
            file_struct['std_cols'].append(c)

    file_struct['n_per_sample']=0
    file_struct['n_BAF']= file_struct['header'].count('B Allele Freq')
    file_struct['n_LRR']= file_struct['header'].count('Log R Ratio')
    file_struct['n_GT'] = file_struct['header'].count('GType')
    assert file_struct['n_BAF'] >= 1, 'No BAF data present'
    file_struct['n_per_sample'] += 1
    assert file_struct['n_LRR'] >= 1, 'No LRR data present'
    file_struct['n_per_sample'] += 1
    assert file_struct['n_BAF'] == file_struct['n_LRR'], "n LRR and n BAF mismatch"
    assert file_struct['n_GT'] == 0 or file_struct['n_GT'] == file_struct['n_BAF'], 'n GType and n BAF mismatch'
    if(file_struct['n_GT']) >=1:
         file_struct['n_per_sample'] += 1
    file_struct['exp_n'] = file_struct['n_std'] + file_struct['n_BAF'] + file_struct['n_LRR'] + file_struct['n_GT']
    file_struct['col_list'] = file_struct['header'].split(',')
    obj.file_struct = file_struct
    

def get_file_order(obj):
    file_order = {}
    file_order['per_sample_cols'] = ['B Allele Freq', 'Log R Ratio', 'GType']
    if obj.file_struct['n_per_sample']==2:
        file_order['per_sample_cols'] = file_order['per_sample_cols'][:-1]
    
    sample_tup_list = list(zip(*[iter(obj.file_struct['col_list'][obj.file_struct['n_std']:])]*obj.file_struct['n_per_sample'])) # zip(*[iter(L)]*2) thankyou stackoverflow: https://stackoverflow.com/questions/23286254/how-to-convert-a-list-to-a-list-of-tuples*2))
    assert len(sample_tup_list) == obj.file_struct['n_BAF']
    # group()
    first_sample = sample_tup_list[0]
    match = re.search(file_order['per_sample_cols'][0], first_sample[0])
    if not match:
        file_order['per_sample_cols'][0], file_order['per_sample_cols'][1] = file_order['per_sample_cols'][1], file_order['per_sample_cols'][0]
    file_order['samples'] = []
    for sample in sample_tup_list:
        assert sample[0].split('.'+file_order['per_sample_cols'][0])[0] == sample[1].split('.'+file_order['per_sample_cols'][1])[0], 'Not all samples follow the same ordering of '+file_order['per_sample_cols'][0]+' '+file_order['per_sample_cols'][1]
        file_order['samples'].append(sample[0].split('.'+file_order['per_sample_cols'][0])[0])
    
    obj.file_order = file_order
    return()


def dedup_samples(obj):
    if np.size(np.unique(obj.file_order['samples'])) == np.size(obj.file_order['samples']):
        return()
    samples = []
    for sample in obj.file_order['samples']:
        if not sample in samples:
            samples.append(sample)
        else:
            n=0
            new_sample = sample
            while new_sample in samples:
                n += 1
                new_sample = sample+':'+str(n)
            samples.append(new_sample)
    obj.file_order['samples'] = samples
    return()


def pl_header(obj):
    type_d = {
        'Name': pl.Utf8,
        'Chr': pl.Categorical,
        'Position': pl.UInt32,
        'GType': pl.Categorical,
        'B Allele Freq': pl.Float32,
        'Log R Ratio': pl.Float32
    }

    col_types = []
    header_cols = []
    for c in obj.file_struct['std_cols']:
        col_types.append(type_d[c])
        header_cols.append(c)
    for s in obj.file_order['samples']:
        for i in range(len(obj.file_order['per_sample_cols'])):
            col_types.append(type_d[obj.file_order['per_sample_cols'][i]])
            header_cols.append(s+'.'+obj.file_order['per_sample_cols'][i])
    obj.col_types = col_types
    obj.header_cols = header_cols
    return()





#----- class defs

# '''class for data to split into n ind chunks'''
class sampleDataPartition():
    def __init__(self, input: str, target_n: int):
        self.input = input
        self.target_n = target_n
        self.df = []
        self.file_struct={}
        self._file_struct()
        self.clean_cols = []
        self._make_clean_cols()
        self.load_data()
        self.samples= []
        self.get_samples()
        self.n_samples = len(self.samples) # total number in input file
        self.n_per_partition=[]
        self.define_partition_n()
        self.prefix = os.path.splitext(input)[0]
    
    def _file_struct(self):
        file_struct(self.input)

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
        self.file_struct={}
        self._file_struct()
        self.clean_cols = []
        self._make_clean_cols()
        self.df = []
        self.load_data()
        self.samples= []
        self.get_samples()
    
    def _file_struct(self.input):
        file_struct(self)

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
    def __init__(self, input: str, geno: float, map_file: None | str):
        self.input = input
        self.geno = geno
        self._make_file_struct()
        self._get_file_order()
        self._dedup_samples()
        self._pl_header()
        self.get_pfb()

    def _make_file_struct(self):
        make_file_struct(self)

    def _get_file_order(self):
        get_file_order(self)

    def _dedup_samples(self):
        dedup_samples(self)

    def _pl_header(self):
        pl_header(self)

    def get_pfb(self):
        col_dict = dict(zip(self.header_cols, self.col_types))
        q = (
            pl.scan_csv(self.input, separator='\t', has_header=False, skip_rows=1, with_column_names=lambda cols: self.header_cols, dtypes = col_dict)
            .select( [*[pl.col(c) for c in  t.file_struct['std_cols']],pl.col("^*.B Allele Freq$") ])
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
            .select([*t.file_struct['std_cols'],"BAF","n_miss","n_call"])
            )
        s= q.collect( streaming=True  )
        s=s.filter(pl.col("n_miss")/(pl.col("n_miss")+pl.col("n_call")) < self.geno)
        s=s.select([
            "Name",
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
