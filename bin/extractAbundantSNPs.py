#!/usr/bin/env python3
"""
filter snps based on threshold of occurence in sequences
"""

import os
import sys
import logging
import argparse
import textwrap
# import pycountry
from textwrap import dedent

import pandas as pd

from utils import mkdir
from utils import fasta_iterator

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent(__doc__)
    )

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--vcf', required=True, type=str,
                                dest="vcf", metavar="<str>",
                                help="path to the VCF file with samples having SNP represented by 1 and that with no "
                                     "snp represented as 0 "
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    required_group.add_argument('--min-freq', type=float,
                                dest="min_freq", metavar="<float>", default=0.2,
                                help="proportion of the minimum number of sequences which contain the SNP"
                                )

    required_group.add_argument('--max-freq', type=float,
                                dest="max_freq", metavar="<float>", default=0.8,
                                help="proportion of the maximum number of sequences which contain the SNP"
                                )
    required_group.add_argument('--snp-type', required=True, type=str,
                                dest="snp_type", metavar="<str>", choices=['singleton', 'multiple', 'conserved', 'all'],
                                default="singleton",
                                help="snp type, singleton (only a single snp per position), multiple (more than one "
                                     "snp per position) and conserved (snps that occur commonly across all the "
                                     "sequences)",
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def extract_most_abundant_snps(vcf, prefix, min_freq, max_freq, snp_type, outdir):
    """

    :param vcf:
    :param accession:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    # read snps from vcf file
    df = pd.read_csv(vcf, sep="\t", skiprows=3)

    # columns to remove when performing SNPs count
    to_remove = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

    # get list of single or multiple snps
    multiple = list()
    snps_dict = {row[1]: row[4] for row in df.values}
    for k, v in snps_dict.items():
        if len(list(v.split(","))) > 1:
            multiple.append(k)

    to_sum = [elem for elem in df.columns if elem not in to_remove]
    min_seqs = int(min_freq * len(to_sum))
    max_seqs = int(max_freq * len(to_sum))

    logging.info("Total number of sequence: {}".format(len(to_sum)))
    logging.info("Minimum number of sequences with {} SNPs: {}".format(snp_type, min_seqs))
    logging.info("Maximum number of sequences with {} SNPs: {}".format(snp_type, max_seqs))

    # separate single snp (singletons) from multiple snps per position
    # only singletons
    if snp_type == 'singleton':
        singleton_df = df[~df['POS'].isin(multiple)]
        # first = singleton_df.columns.get_loc("FORMAT") + 1
        # last = len(list(singleton_df.columns))
        # first_colname = singleton_df.columns[first]
        # last_colname = singleton_df.columns[last-1]
        # singleton_df['seqs'] = singleton_df.loc[:, first_colname:].sum(axis=1)
        # singleton_df['seqs'] = singleton_df.iloc[:, first:].sum(axis=1)
        singleton_df['seqs'] = singleton_df[to_sum].sum(axis=1)
        singleton_df = singleton_df[(singleton_df['seqs'] > min_seqs)].sort_values('seqs', ascending=False)
        logging.info("{} SNPs in {} category".format(len(singleton_df.POS), snp_type))
        result = ' '.join(str(item) for item in sorted(list(singleton_df.POS)))

    # highly variable
    elif snp_type == 'multiple':
        multiple_df = df.loc[df['POS'].isin(multiple)]
        multiple_df['seqs'] = multiple_df[to_sum].sum(axis=1)
        logging.info("{} SNPs in {} category".format(len(multiple_df.POS), snp_type))
        result = ' '.join(str(item) for item in sorted(list(multiple_df.POS)))

    # count snps in highly conserved regions
    elif snp_type == 'conserved':
        conserved_df = df[~df['POS'].isin(multiple)]
        # to_sum = [elem for elem in conserved_df.columns if elem not in to_remove]
        conserved_df['seqs'] = conserved_df[to_sum].sum(axis=1)
        conserved_df = conserved_df[(conserved_df['seqs'] >= max_seqs)].sort_values('seqs', ascending=False)
        logging.info("{} SNPs in {} category".format(len(conserved_df.POS), snp_type))
        result = ' '.join(str(item) for item in sorted(list(conserved_df.POS)))

    # count snps in all the sequences 
    else:
        df['seqs'] = df[to_sum].sum(axis=1)
        df1 = df[(df['seqs'] > min_freq)].sort_values('seqs', ascending=False)
        logging.info("{} SNPs in {} category".format(len(df1.POS), snp_type))
        result = ' '.join(str(item) for item in sorted(list(df1.POS)))

    out_file = os.path.join(outdir, prefix) + '.txt'

    with open(out_file, 'w') as f_obj:
        f_obj.write(result)
    return result


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.vcf is not None and isinstance(args.vcf, str):
        args.vcf = args.vcf
        logging.info("input vcf = {:^10}".format(os.path.basename(args.vcf), str))

    logging.info("extracting snps category: {}".format(args.snp_type))
    extract_most_abundant_snps(vcf=args.vcf, prefix=args.prefix, min_freq=args.min_freq, max_freq=args.max_freq, snp_type=args.snp_type, outdir=args.outdir)


if __name__ == '__main__':
    main()
