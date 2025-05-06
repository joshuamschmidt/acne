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

parser.add_argument('tool', metavar='TOOL', type=str, nargs=1, choices={
                    "partition", "split", "pfb"}, help='which tool to run')

required.add_argument('--input', type=str, dest='input',
                      help='input file with at minimum SNP Name and LRR and BAF cols per sample')

optional.add_argument('--output', type=str,
                      dest='output',
                      help='output:\nif tool is "pfb", then creates a pfb file\nit has no effect for other tools')

optional.add_argument('--n', type=int, dest='n_per_partition',
                      help='how many inds should large input be split into',
                      default=None)

optional.add_argument('--prefix', type=str, dest='prefix',
                      help='batch prefix to attach to ind sample files when splitting',
                      default=None)

optional.add_argument('--geno', type=float, dest='geno',
                      help='SNPs/markers with missingness fraction more than this are excluded from PFB',
                      default=0.02)

optional.add_argument('--samplefilter', type=str,
                      dest='samplefilter',
                      help='samplefilter is 2 two col file listing sample name and binary flag (1 == keep, 0 = filter)',
                      default=None)


'''
TODO: infer geno missingness on GType not BAF
Add ability to set median LRR
also filter samples with optional sample list.
'''
# sub class defs


class fileStructure():
    def __init__(self, file: str):
        self.file = file
        self.__get_header()
        self.__get_std()    
        self.n_BAF = self.__n_data_cols('B Allele Freq')
        self.n_LRR = self.__n_data_cols('Log R Ratio')
        self.n_per_sample = 2
        self.n_GT = self.__n_data_cols('GType')
        if self.n_GT >= 1:
            self.n_per_sample = 3
        self.n_expected = len(self.std_cols) + self.n_per_sample * self.n_BAF
        self.__validate()
        self.size = int(round(os.stat(self.file).st_size / 1e9, 0))

    def __get_header(self):
        with open(self.file, 'rt') as fh:
            self.header = fh.readline().strip().split('\t')
            assert len(self.header) > 1, self.file + ' not TSV?'

    def __get_std(self):
        assert 'Name' in self.header, 'BAF file must have SNP Name column'
        possible_cols = ['Name', 'Chr', 'Position']
        self.std_cols = []
        for p in possible_cols:
            if p in self.header[:3]:
                self.std_cols.append(p)
        assert(len(self.std_cols)) >= 1, 'Error in STD cols'

    def __n_data_cols(self, col):
        return(sum(h.count(col) for h in self.header))

    def __validate(self):
        assert self.n_BAF >= 1, 'Error: No BAF data present'
        assert self.n_LRR >= 1, 'Error: No LRR data present'
        assert self.n_BAF == self.n_LRR, 'Err: n BAF='+self.n_BAF+' n LRR='+self.n_LRR
        assert self.n_GT == 0 or self.n_GT == self.n_BAF, 'n GType and n BAF mismatch'
        assert len(self.header) == self.n_expected, 'more columns than expected'


class sampleOrder():
    def __init__(self, fileStructure, samplefilter: str):
        self.samplefilter = samplefilter
        self.__get_samples(fileStructure)
        self.__dedup_samples()
        self.__filter_samples()
        self.__validate()

    def __get_samples(self, fileStructure):
        sample_tups = list(zip(*[iter(fileStructure.header[len(fileStructure.std_cols):])] * fileStructure.n_per_sample))
        '''
        fold the header into a list of tuples
        size n=n_per_sample.
        idiom: zip(*[iter(L)]*2):
        https://stackoverflow.com/questions/23286254/how-to-convert-a-list-to-a-list-of-tuples*2))
        '''
        first_sample = sample_tups[0]
        self.per_sample_cols = []
        for i in range(fileStructure.n_per_sample):
            self.per_sample_cols.append(first_sample[i].split('.')[1])

        self.samples = []
        for i in range(len(sample_tups)):
            for j in range(len(self.per_sample_cols)):
                if j == 0:
                    this_sample = sample_tups[i][j].split('.'+self.per_sample_cols[j])[0]
                    self.samples.append(this_sample)
                if j >= 1:
                    assert sample_tups[i][j].split('.'+self.per_sample_cols[j])[0] == this_sample, 'Err invalid ordering of sample data: from sample: '+this_sample

    def __dedup_samples(self):
        if np.size(np.unique(self.samples)) == np.size(self.samples):
            self.unique_samples = self.samples
            return()
        self.unique_samples = []
        for sample in self.samples:
            if sample not in self.unique_samples:
                self.unique_samples.append(sample)
            else:
                n = 0
                new_sample = sample
                while new_sample in self.unique_samples:
                    n += 1
                    new_sample = sample+':'+str(n)
                self.unique_samples.append(new_sample)

    def __filter_samples(self):
        self.to_filter = []
        if self.samplefilter is not None:
            with open(self.samplefilter, 'rt') as sh:
                for line in sh:
                    if not line.startswith("Sample_ID"):
                        line = line.strip().split('\t')
                        if int(line[1]) == 0:
                            self.to_filter.append(line[0])
            self.filter_idx=[self.samples.index(f) for f in self.to_filter if f in self.samples]
            self.pl_filter=[self.unique_samples[idx]+'.'+c for c in self.per_sample_cols for idx in self.filter_idx]
        else:
            self.filter_idx = []
            self.pl_filter = []

    def __validate(self):
        assert len(self.samples) == len(self.unique_samples), 'Error when dedup samples'


class plSchema():
    def __init__(self, fileStructure, sampleOrder):
        self.__make_schema(fileStructure, sampleOrder)

    def __make_schema(self, fileStructure, sampleOrder):
        types = {
            'Name': pl.Utf8,
            'Chr': pl.Categorical,
            'Position': pl.UInt32,
            'GType': pl.Categorical,
            'B Allele Freq': pl.Float32,
            'Log R Ratio': pl.Float32
            }
        self.schema = {}

        for c in fileStructure.std_cols:
            self.schema[c] = types[c]

        for s in sampleOrder.unique_samples:
            for c in sampleOrder.per_sample_cols:
                self.schema[s + '.' + c] = types[c]

# ----- class defs


# '''class for data to split into n ind chunks'''
class sampleDataPartition():
    def __init__(self, input: str, prefix: str, target_n: int, samplefilter: str):
        self.input = input
        self.target_n = target_n
        self.samplefilter = samplefilter
        self.fileStructure = fileStructure(self.input)
        self.sampleOrder = sampleOrder(self.fileStructure, self.samplefilter)
        self.plSchema = plSchema(self.fileStructure, self.sampleOrder)
        self.__load_data()
        self.__define_partition_n()
        self.prefix = prefix

    def __load_data(self):
        low_mem = True if self.fileStructure.size > 16 else False

        q = (pl.scan_csv(
            self.input, separator='\t', has_header=False, skip_rows=1,
                        schema=self.plSchema.schema, low_memory = low_mem,
            )
            .drop(self.sampleOrder.pl_filter)
        )
        self.df = q.collect()

    def __define_partition_n(self):
        n_samples = len(self.sampleOrder.unique_samples) - len(self.sampleOrder.filter_idx)
        if(n_samples < 2 * self.target_n):
            self.target_n, remainder = divmod(n_samples, 2)
            if(remainder > 0):
                self.target_n += 1
            print(
                "target_n is larger than input sample_n / 2. Defaulting to splitting input in half")
        n_partitions, remainder = divmod(n_samples, self.target_n)
        self.partition_ns = [self.target_n] * n_partitions
        if(remainder > n_partitions):
            add_all, remainder = divmod(remainder, n_partitions)
            self.partition_ns = [n + add_all for n in self.partition_ns]
        if(remainder > 0):
            for i in range(remainder):
                self.partition_ns[i] += 1

    def write_partition_data(self):
        all_samples = [s for i, s in enumerate(self.sampleOrder.unique_samples) if i not in self.sampleOrder.filter_idx]
        for i, n in enumerate(self.partition_ns):
            samples = all_samples[:n]
            sample_cols = []
            for sample in samples:
                sample_cols += [sample + '.' + col for col in self.sampleOrder.per_sample_cols]
            sub = self.df.select([*self.fileStructure.std_cols, *sample_cols])
            sub.write_csv(self.prefix + "-" + str(i+1) +
                          '.partition', separator='\t')
            all_samples = all_samples[n:]


# '''class for data to split by ind'''
class sampleDataSplit():
    def __init__(self, input: str, prefix: str, samplefilter: str):
        self.input = input
        self.prefix = prefix
        self.samplefilter = samplefilter
        self.fileStructure = fileStructure(self.input)
        self.sampleOrder = sampleOrder(self.fileStructure, self.samplefilter)
        self.plSchema = plSchema(self.fileStructure, self.sampleOrder)
        self.__load_data()

    def __load_data(self):
        q = (
            pl.scan_csv(
                self.input,
                separator='\t',
                has_header=False,
                skip_rows=1,
                with_column_names=lambda cols: list(self.plSchema.schema.keys()), dtypes=self.plSchema.schema,
                )
            .drop(self.sampleOrder.pl_filter)
            )
        self.df = q.collect()

    def write_sample_data(self):
        for i, sample in enumerate(self.sampleOrder.unique_samples):
            if i not in self.sampleOrder.filter_idx:
                sample_cols = [sample + '.' + col for col in self.sampleOrder.per_sample_cols if col != 'GType']
                sub = self.df.select(
                    ["Name", *sample_cols],
                    )
                sub.write_csv(
                    self.prefix + '_' + sample.replace(":", "_") + '.txt',
                    separator='\t',
                    )


# '''class for GtLogRBaf to pfb'''
class pfbObj():
    def __init__(self, input: str, geno: float, samplefilter: str):
        self.input = input
        self.geno = geno
        self.samplefilter = samplefilter
        self.fileStructure = fileStructure(self.input)
        self.sampleOrder = sampleOrder(self.fileStructure, self.samplefilter)
        self.plSchema = plSchema(self.fileStructure, self.sampleOrder)
        self.__get_pfb()

    def __get_pfb(self):
        low_mem = True if self.fileStructure.size > 16 else False
        
        q = (
            pl.scan_csv(self.input, separator='\t', has_header=False, skip_rows=1,
                        schema=self.plSchema.schema, low_memory = low_mem)
            .drop(self.sampleOrder.pl_filter)
            .drop(pl.col("^*.GType$"))
            .select([*[pl.col(c) for c in self.fileStructure.std_cols], pl.col("^*.B Allele Freq$")])
            .with_columns([
                pl.sum_horizontal(pl.col("^*.B Allele Freq$").is_nan()).alias('n_miss'),
                ])
            .with_columns([
                (self.fileStructure.n_BAF - pl.col("n_miss")).alias("n_call")
                ])
            .fill_nan(0)
            .with_columns([
                pl.fold(acc=pl.lit(0),
                        function=lambda acc, x: acc + x, exprs=pl.col("^*.B Allele Freq$")).alias("sum")
            ])
            .with_columns([
                ((pl.col("sum") / (pl.col("n_call") - pl.col("n_miss")))).alias("mean")
            ])
            .fill_nan(0)
            .with_columns([
                ((pl.col("mean")*1000+0.5).cast(pl.Int64)/1000).alias("BAF")
            ])
            .select([*self.fileStructure.std_cols, "BAF", "n_miss", "n_call"])
        )
        s = q.collect()
        s = s.filter(pl.col("n_miss")/( pl.col("n_miss") + pl.col("n_call") ) < self.geno)
        s = s.select([
            *self.fileStructure.std_cols,
            pl.when(pl.col("Name").str.contains("cnv|CNV")).then(
                pl.lit(2)).otherwise(pl.col("BAF")).alias("PFB")
        ])
        self.pfb = s


def main():
    args = parser.parse_args()
    tool = args.tool[0]
    if(tool == 'partition'):
        data = sampleDataPartition(args.input, args.prefix, args.n_per_partition, args.samplefilter)
        data.write_partition_data()

    if(tool == 'split'):
        data = sampleDataSplit(args.input, args.prefix, args.samplefilter)
        data.write_sample_data()

    if(tool == 'pfb'):
        if not args.output:
            parser.error('pfb tool selected: --output must be specified')
        pfb = pfbObj(args.input, args.geno, args.samplefilter)
        pfb.pfb.write_csv(args.output, separator='\t')


if __name__ == '__main__':
    main()
