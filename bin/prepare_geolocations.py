#!/usr/bin/env python3
"""
prepare geo locations input file by merging tab-delimited filw with accessions/taxa and csv file with locations

usage:
python prepare_geolocations.py --metadata <path to tab-delimited text file> --locations <path to locations .csv file> --prefix <M-global> --outDir .
"""

import os
import sys
import logging
import argparse
import textwrap
from textwrap import dedent

import pandas as pd

from utils import mkdir
from utils import str2bool
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
    required_group.add_argument('--metadata', required=True, type=str,
                                dest="metadata", metavar="<str>",
                                help="path to the tab-delimited text file having the column taxa"
                                )
    required_group.add_argument('--locations', required=True, type=str,
                                dest="locations", metavar="<str>",
                                help="full path to the locations (csv) file having the columns: accession,organism,"
                                     "length,strain,country,location,host,date,taxa "
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def prep_geo_locations(metadata, locations, prefix, outdir):
    """

    @param metadata:
    @param locations:
    @param prefix:
    @param outdir:
    @return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    outfile_1 = os.path.join(outdir, prefix) + '_geolocations.csv'
    outfile_2 = os.path.join(outdir, prefix) + '_dates.csv'

    metadata_df = pd.read_csv(metadata, sep="\t")
    locations_df = pd.read_csv(locations, sep=",")

    # merge dataframes
    metadata_df['accession'] = metadata_df.taxa.str.split("|").str[0]
    df = pd.merge(metadata_df, locations_df, on=['accession'])

    cols_to_drop = ['country_y', 'host_y', 'taxa_y']
    res = set(cols_to_drop).issubset(df.columns)

    cols = []
    for col in cols_to_drop:
        if col in df.columns:
            cols.append(col)

    cols_dict = dict()
    for col in cols:
        cols_dict[col.rsplit("_")[0] + '_x'] = col.rsplit("_")[0]

    df = df.drop(cols, axis=1).rename(columns=cols_dict)

    # df['date'] = df.date.str.split('-').str[-1] + '-'+ df.date.str.split('-').str[1] + '-'+ df.date.str.split(
    # '-').str[0]

    df.to_csv(outfile_1, sep=",", index=False)
    df2 = df[['taxa', 'date', 'country', 'location']]
    df2 = df2.rename(columns={'taxa':'accession'})
    df2.to_csv(outfile_2, sep=",", index=False)
    logging.info("output files:\ncsv -> {}".format(outfile_1))
    return outfile_1, outfile_2


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.metadata is not None and isinstance(args.metadata, str):
        logging.info("input metadata = {:^10}".format(args.metadata, str))
    if args.locations is not None and isinstance(args.locations, str):
        logging.info("input locations = {:^10}".format(args.locations, str))

    logging.info("preparing geo locations")
    prep_geo_locations(metadata=args.metadata, locations=args.locations, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
