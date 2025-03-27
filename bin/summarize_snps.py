#!/usr/bin/env python3
"""
summarize SNPs per grouping
"""

import os
import re
import sys
import csv
import logging
import argparse
import itertools
# from itertools import izip_longest
# from itertools import chain
from textwrap import dedent
from collections import Counter
import pandas as pd
import numpy as np

# import pycountry
# from Bio import SeqIO
# from Bio import Entrez

from utils import mkdir
from utils import str2bool

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent('''summarize SNPs per grouping''')
    )
    # helpstr = """python3 fetch_segments.py [options]"""

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--snps', required=True, type=str,
                                dest="snps", metavar="<str>",
                                help="comma separated values (.csv) file containing the columns 'record', 'snps', "
                                     "and 'num_snps' from snipit "
                                )
    required_group.add_argument('--metadata', required=True, type=str,
                                dest="metadata", metavar="<str>",
                                help="path to the file containing metadata"
                                )
    parser.add_argument('--recomb', metavar="<str>", nargs="+", dest="recomb", default='', required=False,
                        help="One or more accessions comma-separated list of recombinant sequences. eg "
                             "NC_014398,NC_014399")
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outdir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def summarise_snps(snps, metadata, recomb, prefix, outdir):
    """

    :param snps:
    :param metadata:
    :param recomb
    :return:
    """

    logging.info("Recombinants: {}".format(recomb))

    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out1 = os.path.join(outdir, prefix) + '.all.csv'
    out2 = os.path.join(outdir, prefix) + '.per.lineage.csv'

    snps_df = pd.read_csv(snps, sep=",")
    metadata_df = pd.read_csv(metadata, sep="\t")

    df2 = pd.merge(left=snps_df, right=metadata_df, left_on='record', right_on='accession')

    # print(df2[(df2.strain_type == 'vaccine')][['accession', 'Lineage', 'strain', 'strain_type']])

    df3 = df2[['accession', 'lineage', 'date', 'strain_type', 'snps']]
    if len(recomb) > 0:
        df3 = df3[~df3['accession'].isin(recomb)]

    df3['positions'] = df3['snps'].apply(lambda x: re.findall(r"[-+]?\d*\.\d+|\d+", str(x))).map(
        lambda x: [int(i) if '.' not in i else float(i) for i in x]).tolist()
    df3['SNPs'] = df3['snps'].apply(lambda x: re.findall(r"[ATCG]{2}", str(x))).map(
        lambda x: [str(i) if '.' not in i else str(i) for i in x]).tolist()

    # group by lineage and strain type (vaccine or non-vaccine)
    grps = df3.groupby(['lineage', 'strain_type'])
    d = dict()
    for name, grp in grps:
        snps_lists = grp['SNPs'].tolist()
        snps_list = list(itertools.chain(*snps_lists))
        d[name[0] + ' ' + name[1]] = snps_list

    df4 = pd.DataFrame(d.items(), columns=['name', 'SNPs'])
    df4[['Lineage', 'Strain']] = df4.name.str.split(" ", expand=True)
    df4 = df4[['Lineage', 'Strain', 'SNPs']]
    # print(df4)

    d2 = dict()
    for index, row in df4.iterrows():
        lineage = row['Lineage']
        strain = row['Strain']
        snps = sorted(row['SNPs'])
        mutations = list(Counter(snps).keys())
        counts = list(Counter(snps).values())
        d1 = dict(zip(mutations, counts))
        d2[lineage + ':' + strain] = d1
    df5 = pd.DataFrame(d2)
    # print(df5)
    df6 = df5.T.reset_index().rename(columns={'index': 'SNP'})

    ti = ['AG', 'GA', 'CT', 'TC']
    ts = ['AC', 'AT', 'CA', 'CG', 'GC', 'GT', 'TA', 'TG']
    mutations_cols = list(df6.columns)[1:]

    tv_cols = sorted(list(set(mutations_cols) - set(ti)))
    ti_cols = sorted(list(set(mutations_cols) - set(ts)))

    df6['Ti'] = df6[ti_cols].sum(axis=1)
    df6['Tv'] = df6[tv_cols].sum(axis=1)
    df6['Ti/Tv'] = df6['Ti'].div(df6['Tv'].values)
    # df6.loc[~np.isfinite(df6['Ti/Tv']), 'result'] = np.nan
    df6['Total'] = df6[['Ti', 'Tv']].sum(axis=1)
    df6 = df6.round(1)
    # print(df6)
    df6.to_csv(out2, index=False, header=True)

    # all or individual accessions
    d2_dict = dict()
    for index, row in df3.iterrows():
        accession = row['accession']
        lineage = row['lineage']
        strain = row['strain_type']
        snps = sorted(row['SNPs'])
        mutations = list(Counter(snps).keys())
        counts = list(Counter(snps).values())
        d2 = dict(zip(mutations, counts))
        d2_dict[str(accession) + ':' + lineage + ':' + strain] = d2
    df7 = pd.DataFrame(d2_dict)
    df8 = df7.T.reset_index().rename(columns={'index': 'SNP'})
    df8['Ti'] = df8[ti_cols].sum(axis=1)
    df8['Tv'] = df8[tv_cols].sum(axis=1)
    df8['Ti/Tv'] = df8['Ti'].div(df8['Tv'].values)
    df8['Total'] = df8[['Ti', 'Tv']].sum(axis=1)
    df8 = df8.round(1)
    df8.to_csv(out1, index=False, header=True)


def main():
    parser = parse_args()
    args = parser.parse_args()

    summarise_snps(snps=args.snps, metadata=args.metadata, recomb=args.recomb,
                   prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
