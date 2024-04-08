#!/usr/bin/env python
import argparse
import textwrap

'''required and optional argument parser'''

parser = argparse.ArgumentParser(prog='sample_sheet_check',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                        Checks the format of the input sample_sheet.
                        Adds paths to dummy files if not passed as variable
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


required.add_argument('--samplesheet', type=str, dest='samplesheet',
                      help='input samplesheet')

required.add_argument('--outfile', type=str,
                      dest='outfile',
                      help='validated/updated sample sheet')


def main():
    args = parser.parse_args()
    req_header_cols = ['id', 'gsfile', 'gs_sample_sheet', 'sample_include', 'snp_include']
    req_body = [True, True, True, False, False]

    fout = open(args.outfile, 'wt')

    with open(args.samplesheet, 'rt') as fin:
        req_col_dict = {}
        for i, line in enumerate(fin):
            if i == 0:
                line = line.strip().split(',')
                for col in req_header_cols:
                    assert col in line, 'Header col ' + col + ' is missing!'
                    req_col_dict[col] = line.index(col)
                print(','.join(req_header_cols), file=fout)
            else:
                line = line.strip().split(',')
                outline = []
                for i, col in enumerate(req_header_cols):
                    col_val = line[req_col_dict[col]]
                    if(req_body[i]):
                        assert col_val != '', col + ' is an empty string!'
                        outline.append(col_val)
                    else:
                        if col_val == '':
                            if col == 'sample_include':
                                outline.append('./assets/sample_include_dummy.txt')
                            elif col == 'snp_include':
                                outline.append('./assets/snp_include_dummy.txt')
                        else:
                            outline.append(col_val)
                print(','.join(outline), file=fout)
    fout.close()


if __name__ == '__main__':
    main()
