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
                      help='output:\nif tool is "pfb", then creates a pfb file\nit has no effect for other tools')

optional.add_argument('--samplesheet', type=str, dest='samples',
                      help='genome studio sample sheet. required to update sample names from sentrix barcode + pos format',
                      default=None)

optional.add_argument('--genorate', type=float, dest='geno',
                      help='batch prefix to attach to ind sample files when splitting',
                      default=None)

optional.add_argument('--lrrsd', type=float, dest='lrrsd',
                      help='LRR deviation threshold (corresponds to qclrrsd in PENN CNV)',
                      default=None)

