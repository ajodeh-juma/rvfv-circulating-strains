#!/usr/bin/env python3
"""
prepare consensus genomes from sequenced samples for phylogenetic analysis

usage:
python prepare_sequenced_samples.py --consensus-dir <fasta-dir> --metadata <metadata> --outDir .
"""

import os
import sys
import logging
import pycountry
import argparse
import textwrap
from textwrap import dedent
from datetime import datetime

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
    required_group.add_argument('--consensus-dir', required=True, type=str,
                                dest="consensus", metavar="<str>",
                                help="path to the directory containing consensus genomes in fasta format"
                                )
    required_group.add_argument('--metadata', required=True, type=str,
                                dest="metadata", metavar="<str>",
                                help="path to the comma-separated values file containing samples metadata"
                                )
    required_group.add_argument('--summary', required=True, type=str,
                                dest="summary", metavar="<str>",
                                help="path to the comma-separated values file containing samples quality on genome "
                                     "coverage "
                                )
    parser.add_argument('--min-cov', required=False, type=float, default=80.0,
                        dest="min_cov", metavar="<float>",
                        help="minimum genome coverage"
                        )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def prep_data(consensus, metadata, summary, prefix, min_cov, outdir):
    """
    search and download complete rift valley fever virus sequences from NCBI

    :param consensus
    :param metadata:
    :@param summary
    :@param: prefix
    :@param: min_cov
    :param outdir:
    :return:
    """

    # use module pycountry to create a dictionary of counties and their 3-letter codes/symbols
    countries = dict()
    t = list(pycountry.countries)
    for country in t:
        countries[country.name.split(',')[0]] = country.alpha_3

    deposited = ['DVS-230', 'DVS-333', 'DVS-356', 'DVS-321', 'DVS-372']

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_file = os.path.join(outdir, prefix) + '.seq.samples.csv'
    out_fasta = os.path.join(outdir, prefix) + '.seq.samples.fasta'

    # read metadata
    metadata_df = pd.read_csv(metadata, sep=",")
    # read consensus genome quality
    summary_df = pd.read_csv(summary, sep=",")
    # merge
    df = pd.merge(summary_df, metadata_df, on="sample_name")

    # filter on coverage
    df = df[(df.pct_covered_bases > float(min_cov))]

    organism = 'Rift Valley fever virus'
    df['organism'] = organism
    df = df[
        ['sample_name', 'organism', 'longest_no_N_run', 'sample_name', 'country', 'location', 'host', 'date', 'fasta']]
    df.columns = ['accession', 'organism', 'length', 'strain', 'country', 'location', 'host', 'date', 'fasta']

    # exclude deposited sequences
    df = df[~df['accession'].isin(deposited)]

    # create a dict
    d = df.set_index('accession').T.to_dict('list')

    df = df[['accession', 'organism', 'length', 'strain', 'country', 'location', 'host', 'date']]
    df['country'] = df['country'].map(countries)
    df['host'] = df['host'].str.lower()
    # print(df.host.unique())
    # logging.info("dataframe from sequenced samples\n {}".format(df))

    df.to_csv(out_file, index=False)

    with open(out_fasta, 'w') as f_obj:
        for k, v in d.items():
            accession = k
            host = v[5].lower()
            date = v[6]
            date = datetime.strptime(date, '%d-%m-%Y').date()
            country = v[3]
            location = v[4]
            country = countries.get(country)

            if accession in deposited:
                # print(accession)
                continue
            else:
                new_accession = '>{}|{}|{}|{}|{}'.format(accession, host, country, location, date)
                fasta = os.path.join(consensus, v[7])

            if os.path.exists(fasta):
                for line in open(fasta):
                    if line.startswith('>'):
                        f_obj.write(new_accession + "\n")
                    else:
                        f_obj.write(textwrap.fill(line, 80))
                        f_obj.write("\n")
    logging.info("{} sequences retrieved with specified min coverage {}".format(len(list(df.accession)), min_cov))
    logging.info("output files:\ncsv -> {}\nfasta -> {}".format(out_file, out_fasta))
    return out_file, out_fasta


def main():
    parser = parse_args()
    args = parser.parse_args()

    prep_data(consensus=args.consensus, metadata=args.metadata, summary=args.summary, prefix=args.prefix,
              min_cov=args.min_cov,
              outdir=args.outdir)


if __name__ == '__main__':
    main()
