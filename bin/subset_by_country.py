#!/usr/bin/env python3
"""
fetch viral sequences from genbank using Biopython
"""

import os
import sys
import logging
import argparse
import textwrap
import pycountry
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
    required_group.add_argument('--fasta', required=True, type=str,
                                dest="fasta", metavar="<str>",
                                help="path to the FASTA format file"
                                )
    required_group.add_argument('--csv', required=True, type=str,
                                dest="csv", metavar="<str>",
                                help="path to the CSV file containing the FASTA sequences metadata"
                                )
    required_group.add_argument('--countries', required=True, type=str, nargs="+",
                                dest="countries", metavar="<str>",
                                help="countries to subset sequences"
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    # required_group.add_argument('--out-fasta', metavar="<file>", type=str, dest="out_fasta", help="output file to write sequences in fasta format")
    # required_group.add_argument('--out-csv', metavar="<file>", type=str, dest="out_csv", help="output file to write metadata in csv format")
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def subset_data(fasta, csv, countries, prefix, outdir):
    """
    search and download complete rift valley fever virus sequences from NCBI

    :param fasta:
    :param csv:
    :param countries:
    :param prefix:
    :param out_fasta:
    :param out_csv:
    :param outdir:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_csv = os.path.join(outdir, prefix) + '.csv'
    out_fasta = os.path.join(outdir, prefix) + '.fasta'

    df = pd.read_csv(csv, sep=",")

    # use module pycountry to create a dictionary of counties and their 3-letter codes/symbols
    countries_d = dict()
    t = list(pycountry.countries)
    for c in t:
        countries_d[c.name.split(',')[0]] = c.alpha_3
    codes = [countries_d.get(country) for country in countries]

    # subset countries
    country_df = df[df['country'].isin(codes)]
    country_df.to_csv(out_csv, sep=",", index=False)

    accessions = list(country_df.accession)
    logging.info("found {} sequences".format(len(accessions)))

    # iterate through the input fasta file
    fasta_dict = dict(fasta_iterator(fasta))

    f = []

    for k, v in fasta_dict.items():
        if k.split('|')[0] in accessions:
            f.append(k)

    with open(out_fasta, 'w') as outfile:
        for acc in f:
            outfile.write(">" + acc + "\n")
            outfile.write(textwrap.fill(fasta_dict.get(acc), 80))
            outfile.write("\n")
    return out_csv, out_fasta


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.fasta is not None and isinstance(args.fasta, str):
        args.fasta = args.fasta
        logging.info("input fasta = {:^10}".format(args.fasta, str))

    if args.csv is not None and isinstance(args.csv, str):
        args.csv = args.csv
        logging.info("input csv = {:^10}".format(args.csv, str))

    if args.countries is not None:
        countries = ' '.join(str(k) for k in args.countries)
        logging.info("selected countries = {:^10}".format(countries))

    subset_data(fasta=args.fasta, csv=args.csv, countries=args.countries, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
