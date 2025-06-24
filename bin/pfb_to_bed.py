#!/usr/bin/env python
import argparse
import textwrap
import numpy as np

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='pfb_to_bed',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                        create a bed format file from a penn cnv pfb file
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


required.add_argument('--pfb', type=str, dest='pfb',
                      help='pfb file with name, chr, pos')

optional.add_argument('--window', type=int, dest='window',
                      help='window size in bp',
                      default=500000)

required.add_argument('--chr_sizes', type=str, dest='chr_sizes',
                      help='list of chr sizes form ucsc')

required.add_argument('--outfile', type=str,
                      dest='outfile',
                      help='name of output file')


def main():
    args = parser.parse_args()
    fout = open(args.outfile, 'wt')

    chr_sizes = {}
    ucsc_prefix = False
    with open(args.chr_sizes, 'rt') as cin:
        for line in cin:
            line = line.strip().split('\t')
            chr_sizes[str(line[0])] = int(line[1])
    if all([k.startswith('chr') for k in chr_sizes.keys()]):
        ucsc_prefix = True
    with open(args.pfb, 'rt') as fin:
        pfb_prefix = False
        line_info = {}
        for i, line in enumerate(fin):
            if i == 0:
                line = line.strip().split('\t')
                for j, l in enumerate(line):
                    line_info[l] = j
            else:
                line = line.strip().split('\t')
                chrom = line[line_info['Chr']]
                if i == 1:
                    if chrom.startswith('chr'):
                        pfb_prefix = True
                    if ucsc_prefix and pfb_prefix != ucsc_prefix:
                        chr_sizes = {key.replace("chr", ""): value for key, value in chr_sizes.items()}
                chrom_size = chr_sizes[str(chrom)]
                name = line[line_info['Name']]
                pos = int(line[line_info['Position']])

                start = np.max([pos - 1 - args.window, 0])
                end = np.min([pos + args.window, chrom_size])
                print(chrom + '\t' + str(start) + '\t' + str(end) + '\t' + name, file=fout)
    fout.close()


if __name__ == '__main__':
    main()
