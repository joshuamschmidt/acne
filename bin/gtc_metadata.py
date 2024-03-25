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

required.add_argument('--outfile', type=str,
                      dest='output',
                      help='name of output file')

optional.add_argument('--samplesheet', type=str, dest='samplesheet',
                      help='genome studio sample sheet. required to update sample names from sentrix barcode + pos format',
                      default=None)

optional.add_argument('--callrate', type=float, dest='callrate',
                      help='minimum callrate to include sample',
                      default=0.90)

optional.add_argument('--lrrsd', type=float, dest='lrrsd',
                      help='LRR deviation threshold (corresponds to qclrrsd in PENN CNV)',
                      default=0.2)


def _not_header_line(s):
    # return true if a line equals "Sample_ID"
    return not s.startswith("Sample_ID")


def main():
    args = parser.parse_args()
    if(args.samplesheet is not None):
        samples = {}
        cols = {}
        with open(args.samplesheet, 'rt') as sh:
            for line in dropwhile(_not_header_line, sh):
                if line.startswith("Sample_ID"):
                    line = line.strip().split(',')
                    cols['barcode'] = line.index('SentrixBarcode_A')
                    cols['pos'] = line.index('SentrixPosition_A')
                    cols['sample'] = line.index('Sample_ID')
                    cols['group'] = line.index('Sample_Group')
                else:
                    line = line.strip().split(',')
                    samples[line[cols['barcode']] + '_' + line[cols['pos']]] = line[cols['group']] + ':' + line[cols['sample']]

    fout = open(args.outfile, 'wt')
    with open(args.metadata, 'rt') as mh:
        cols = {}
        n_line = 0
        fout.write("Sample_ID"+'\t'+"Keep"+'\n')
        for line in mh:
            n_line += 1
            line = line.strip().split('\t')
            if(n_line == 1):
                assert 'gtc' in line, 'Error reading gtc metadata'
                cols['gtc'] = line.index('gtc')
                cols['callrate'] = line.index('call_rate')
                cols['lrrsd'] = line.index('logr_deviation')
            else:
                keep = 1
                gtc = line[cols['gtc']].split('.gtc')[0]
                sample = gtc
                if(args.samplesheet is not None):
                    sample = samples[gtc]
                callrate = float(line[cols['callrate']])
                lrrsd = float(line[cols['lrrsd']])
                if callrate < args.callrate or lrrsd > args.lrrsd:
                    keep = 0
                fout.write(sample+'\t'+str(keep)+'\n')
    fout.close()


if __name__ == '__main__':
    main()
