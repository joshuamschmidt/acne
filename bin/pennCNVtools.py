#!/usr/bin/env python
import argparse
import numpy as np
import polars as pl
import textwrap
import os

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

parser.add_argument('tool',metavar='TOOL', type=str, nargs=1, choices={"partition", "pfb", "split"}, help='which tool to run')

required.add_argument('--input', type=str, dest='input',
                    help='input file from genome studio with GT, LogR and BAF cols per sample')

optional.add_argument('--output', type=str,
                       dest='output',
                       help='output:\nif tool is "pfb", then creates a pfb file\nif "split", prefix for breaking into n chunks of inds\n if "split-ind", it has no effect')

optional.add_argument('--n', type=int, dest='n_per_partition',
                      help='how many inds should large input be split into',
                      default=None)


#----- class defs

# '''class for data to split into n ind chunks'''
class sampleDataParition():
    def __init__(self, input: str, target_n: int):
        self.input = input
        self.target_n = target_n
        self.df = []
        self.load_data()
        self.samples= []
        self.get_samples()
        self.n_samples = len(self.samples) # total number in input file
        self.n_per_partition=[]
        self.define_partition_n()
        self.prefix = os.path.splitext(input)[0]

    def load_data(self):
        q = pl.scan_csv(self.input, sep='\t')
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
            sub.write_csv(self.prefix + '_partition-' + str(j+1) + '.txt', sep='\t')
            all_samples = all_samples[n_part:]




# '''class for data to split by ind'''
class sampleData():
    def __init__(self, input: str,):
        self.input = input
        self.df = []
        self.load_data()
        self.samples= []
        self.get_samples()

    def load_data(self):
        q = pl.scan_csv(self.input, sep='\t')
        self.df = q.collect()
    
    def get_samples(self):
        self.samples=self.df.select(pl.col("^*.GType$")).columns
        self.samples=[s.split('.GType')[0] for s in self.samples]
    
    def write_sample_data(self):
        for s in self.samples:
            col_1, col_2, col_3 = s+".GType", s+".Log R Ratio", s+".B Allele Freq"
            sub=self.df.select(["Name", "Chr", "Position", col_1, col_2, col_3])
            sub.write_csv(s+'.input.txt', sep='\t')


# '''class for GtLogRBaf to pfb'''
class pafObj():
    def __init__(self, input: str,):
        self.input = input
        self.get_paf()
    def get_paf(self):
        q = (
            pl.scan_csv(self.input, sep='\t')
            .select(["Name", "Chr", "Position", pl.col("^*.GType$")])
            .with_columns([
                pl.sum(pl.all().exclude(["Name", "Chr", "Position"]) == "AA").alias('ref'),
                pl.sum(pl.all().exclude(["Name", "Chr", "Position"]) == "AB").alias('het'),
                pl.sum(pl.all().exclude(["Name", "Chr", "Position"]) == "BB").alias('alt'),
                ])
            .with_columns([
                ( (pl.col("het") + (pl.col("alt")*2)) / ((pl.col("ref") + pl.col("het") + pl.col("alt")) * 2)).alias("PFB")
                ])
            .select(["Name", "Chr", "Position","PFB"])
            )
        df = q.collect()
        self.paf = df

def main():
    args = parser.parse_args()
    tool=args.tool[0]
    if(tool=='partition'):
        data = sampleDataParition(args.input,args.n_per_partition)
        data.write_partition_data()

    if(tool=='pfb'):
        paf = pafObj(args.input)
        paf.paf.write_csv(args.output, sep='\t')

    if(tool=='split'):
        data = sampleData(args.input)
        data.write_sample_data()

if __name__ == '__main__':
    main()
