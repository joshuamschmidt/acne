#!/usr/bin/env python
import argparse
import textwrap
import numpy as np

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='extract_gc',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                        created a 2 col file of SNP Name and GC content
                        extracts from out put of bedtools nuc
                        note scale of bedtools nuc is 0-1, but
                        pennCNV requires 0-100
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


required.add_argument('--nuc', type=str, dest='nuc',
                      help='bedtools nuc output')

required.add_argument('--outfile', type=str,
                      dest='outfile',
                      help='name of output file')


def main():
    args = parser.parse_args()
    fout = open(args.outfile, 'wt')

    get_cols = ['4_usercol', '6_pct_gc']

    with open(args.nuc, 'rt') as fin:
        col_dict = {}
        line_info = {}
        for i, line in enumerate(fin):
            if i == 0:
                print('Name' + '\t' + 'GC', file=fout)
                line = line.strip().split('\t')
                for col in get_cols:
                    col_dict[col] = line.index(col)
            else:
                line = line.strip().split('\t')
                snp_name = line[col_dict['4_usercol']]
                snp_gc = float(line[col_dict['6_pct_gc']]) * 100
                snp_gc_format = "{:.8f}".format(snp_gc)
                print(snp_name + '\t' + str(snp_gc_format), file=fout)
    fout.close()


if __name__ == '__main__':
    main()
