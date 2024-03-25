#!/usr/bin/env python
import argparse
import os
import sys
import textwrap
from itertools import dropwhile, takewhile

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='gtc_metadata',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                        parse gtc metadata file created by bcftools +gts2vcf
                        and create file for sample filtering
                        ------------------------------------------------
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


required.add_argument('--metadata', type=str, dest='metadata',
                      help='gtc metadata')

required.add_argument('--output', type=str,
                      dest='output',
                      help='name of output file')

optional.add_argument('--samplesheet', type=str, dest='samplesheet',
                      help='genome studio sample sheet. required to update sample names from sentrix barcode + pos format',
                      default=None)

optional.add_argument('--callrate', type=float, dest='callrate',
                      help='minimum callrate to include sample',
                      default=0.0)

optional.add_argument('--lrrsd', type=float, dest='lrrsd',
                      help='LRR deviation threshold (corresponds to qclrrsd in PENN CNV)',
                      default=None)


def _not_header_line(s):
    # return true if a line equals "Sample_ID"
    return not s.startswith("Sample_ID")


def main():
    args = parser.parse_args()
    if(args.samplesheet is not None):
        samples = {}
        cols = {}
        with open(args.samplesheet, 'rt') as fh:
            for line in dropwhile(_not_header_line, fh):
                if line.startswith("Sample_ID"):
                    line = line.strip().split(',')
                    cols['barcode'] = line.index('SentrixBarcode_A')
                    cols['position'] = line.index('SentrixPosition_A')
                    cols['sample'] = line.index('Sample_ID')
                else:
                    line = line.strip().split(',')
                    samples[line[cols['barcode']] + '_' + line[cols['position']]] = line[cols['sample']]



if __name__ == '__main__':
    main()
