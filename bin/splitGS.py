#!/usr/bin/env python
import argparse
import numpy as np
import polars as pl
import textwrap

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='splitGS',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent('''\
                        split Genome Studio Output to ind files
                        ------------------------------------------------
                        Generates a ind file for PENN CNV

                        Output Column headers are:
                        Name    Chr     Position        Sample.GT sample.`Log R Ratio` sample 
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

required.add_argument('--input', type=str, dest='input',
                    help='input file from genome studio with GT, LogR and BAF cols per sample')


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


def main():
    args = parser.parse_args()
    data = sampleData(args.input)
    data.write_sample_data()

if __name__ == '__main__':
    main()