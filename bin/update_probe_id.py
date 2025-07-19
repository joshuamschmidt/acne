#!/usr/bin/env python
import argparse
import numpy as np
import polars as pl
import textwrap
import os

from pennCNVtools import fileStructure, plSchema, sampleOrder

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='update_probe_names',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                        Python Utility for updating the Name attr of probes
                        in any type of file for pennCNV
                        ------------------------------------------------
                        With a supplied tab delimited key file, old_name new_name,
                        It is an easy 1-1 swap. Without a key file, it replaces Name with
                        Chr + Pos                                                                             


                        Option "pfb" Generates a PFB file for PENN CNV

                        Output Column headers are:
                        Name    Chr     Position        PFB

                        Option "split" Generates a ind file for
                        PENN CNV

                        Output Column headers are:
                        Name  Sample.GT sample.Log R Ratio Sample.BAF

                        Note - `Join's` are inner -> so we drop probes in the input
                        that do not have a 'Name' appearing in the sanme_file.
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
                      help='input file with at minimum SNP Name and LRR and BAF cols per sample')

required.add_argument('--file-type', type=str, dest='file_type',
                      help='Type of input file: "gc_model", "pfb", or "BAF"')

required.add_argument('--output', type=str,
                      dest='output',
                      help='output:\nif tool is "pfb", then creates a pfb file\nit has no effect for other tools')

optional.add_argument('--name-file', type=str, dest='name_file',
                      help='file containing Name -> New_Name mapping. 2 cols',
                      default=None)


class gcModelFile():
    def __init__(self, input: str, name_file: str, output: str):
        self.input = input
        self.name_file = name_file
        self.output = output
        self.input_plSchema = {
            'Name': pl.Utf8,
            'GC': pl.Float64,
        }
        self.name_file_plSchema = {
            'Name': pl.Utf8,
            'New_Name': pl.Utf8,
        }
    
    def process_file(self):
        if self.name_file is None:
             raise ValueError("A name_file must be used to update a GC model file")
        input = (
            pl.scan_csv(self.input, separator='\t', has_header=True,
                        schema=self.input_plSchema)
        )
        if self.name_file is not None:
                name_file = (
                     pl.scan_csv(self.name_file, separator='\t', has_header=True,
                        schema=self.name_file_plSchema)
                )
                input = input.left_join(name_file, on="Name")
                input.update(name_file, left_on="Name", right_on="Old_Name", how="full", include_nulls=True)
                .sink_csv(self.output, separator='\t') 
            )
            del(q)


def main():
    args = parser.parse_args()


if __name__ == '__main__':
    main()
