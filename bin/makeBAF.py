#!/usr/bin/env python
import argparse
import numpy as np
import polars as pl
import textwrap

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='makePFB',
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        description=textwrap.dedent('''\
                        make PFB from Genome Studio Output
                        ------------------------------------------------
                        Generates a PFB file for PENN CNV

                        Output Column headers are:
                        Name    Chr     Position        PFB 
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

required.add_argument('--output', type=str,
                       dest='output',
                       help='output pfb file')

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
    paf = pafObj(args.input)
    paf.paf.write_csv(args.output, sep='\t')

if __name__ == '__main__':
    main()