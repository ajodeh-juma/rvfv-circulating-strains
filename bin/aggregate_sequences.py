#!/usr/bin/env python3
"""
aggregate sequences form NCBI and locally sequenced ones and merge metadata dataframes

usage:
python aggregate_sequences.py --data-dir <data-dir> --prefix <prefix> --outDir .
"""

import os
import sys
import logging
import argparse
import textwrap
import operator
from textwrap import dedent
from collections import OrderedDict

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
    required_group.add_argument('--data-dir', required=True, type=str,
                                dest="datadir", metavar="<str>",
                                help="path to the directory containing multi-fasta sequences and related metadata"
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def aggregate_data(datadir, prefix, outdir):
    """
    aggregate data both sequence multifasta files and metadata files

    :param datadir
    :param metadata:
    :@param summary
    :@param: prefix
    :@param: min_cov
    :param outdir:
    :return:
    """

    # create output directory if not exists
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out_file = os.path.join(outdir, prefix) + '.csv'
    out_fasta = os.path.join(outdir, prefix) + '.fasta'

    d = dict()
    # loop through the directory
    for fn in os.listdir(datadir):
        basename = fn.split('.')[0]
        if fn.endswith('.csv'):
            csv = os.path.join(datadir, fn)
            if basename not in d and os.path.exists(csv):
                d[basename] = [csv]
            else:
                d[basename].append(csv)
        if fn.endswith('.fasta'):
            fasta = os.path.join(datadir, fn)
            if basename not in d and os.path.exists(fasta):
                d[basename] = [fasta]
            else:
                d[basename].append(fasta)

    f_dict = dict()
    dfs = []
    seq_dict = dict()
    # df = pd.DataFrame()
    # to_deposit = []

    print("Run\tSample\tPercentage(N)\tCoverage\n")
    for k, v in d.items():
        # iterate through the input fasta file
        csv, fasta = sorted(v)[0], sorted(v)[1]
        fasta_dict = dict(fasta_iterator(fasta))

        for header, seq in fasta_dict.items():
            # extract positions with N's
            base_positions = [index for index, base in enumerate(seq.lower()) if base != "n"]
            n_positions = [index for index, base in enumerate(seq.lower()) if base == "n"]
            # compute percentage of N's
            perc_n = str(round(len(n_positions) / len(seq) * 100, 2))
            perc_cov = str(round(len(base_positions) / len(seq) * 100, 2))
            accession = str(header.split('|')[0])
            host = header.split('|')[1]
            organism = 'Rift Valley fever virus'
            length = len(seq)
            country = header.split('|')[2]
            location = header.split('|')[3]
            date = header.split('|')[-1]
            f_dict[accession] = [accession, host, length, country, location, date, perc_n]

            columns = ['accession', 'host', 'length', 'country', 'location', 'date', 'perc_n']
            df = pd.DataFrame.from_dict(f_dict, orient='index', columns=columns)
            # create a taxa column in the dataframe
            fasta_headers = ['accession', 'host', 'country', 'location', 'date', 'perc_n']
            df['taxa'] = df.accession.str.cat([df[k] for k in fasta_headers[1:]], sep="|")

            # create a dict with new headers and sequences
            new_header = accession + '|' + host + '|' + country + '|' + location + '|' + date + '|' + perc_n
            seq_dict[new_header] = seq
            to_deposit_ncbi_dict = dict()
            if float(perc_n) > 0.0:
                print(k, '\t', header, '\t', perc_n, '\t', perc_cov)

            # if float(perc_n) > 0.0 and k != prefix:
            #     to_deposit_ncbi_dict[header] = [k, accession, host, organism, length, perc_n, perc_cov, country, location, date]
            #     cols = ['runid', 'accession', 'host', 'organism', 'length', 'perc_n', 'perc_cov', 'country', 'location', 'date']
            #     to_deposit_df = pd.DataFrame.from_dict(to_deposit_ncbi_dict, orient='index', columns=cols)
            #     to_deposit.append(to_deposit_df)

        # read csv metadata file
        csv_df = pd.read_csv(csv, sep=",")
        csv_df = csv_df.astype({"accession": str})
        if len(csv_df) == 0:
            continue
        else:
            dfs.append(csv_df)

    # to_deposit_data = pd.concat(to_deposit)
    # to_deposit_data = to_deposit_data.sort_values(by=['accession', 'length', 'perc_n'])
    # to_deposit_data = to_deposit_data.drop_duplicates(subset=['accession'], keep='last').sort_values(by=['country'], ascending=True)
    # out2 = os.path.join(outdir, prefix) + '_to_deposit_ncbi.csv'
    # to_deposit_data.to_csv(out2, index=False)

    # concatenate the dataframes
    data = pd.concat(dfs)
    # merge two dataframes and filter on percentage N
    merged = pd.merge(df, data, on="accession")
    merged = merged.sort_values(by=['accession', 'length_y', 'perc_n'])
    duplicate = merged[merged.duplicated('accession')]

    # print(duplicate)
    logging.info("duplicated sequences {}".format(','.join(list(duplicate.accession))))
    merged = merged.drop_duplicates(subset=['accession'], keep='last').sort_values(by=['country_x'], ascending=True)

    # accession,organism,length,strain,country,location,host,date
    with open(out_fasta, 'w') as f_obj:
        for taxa in list(merged.taxa):
            if seq_dict.get(taxa) is None:
                continue
            else:
                header, seq = '|'.join(taxa.split('|')[:-1]), seq_dict.get(taxa)
                f_obj.write('>' + header)
                f_obj.write("\n")
                f_obj.write(textwrap.fill(seq, 80))
                f_obj.write("\n")

    merged = merged[['accession', 'organism', 'length_x', 'strain', 'country_x', 'location_y', 'host_x', 'date_x']]
    merged.columns = ['accession', 'organism', 'length', 'strain', 'country', 'location', 'host', 'date']
    merged = merged.sort_values(by=['country'], ascending=True)
    merged.to_csv(out_file, index=False)
    logging.info("output files:\ncsv -> {}\nfasta -> {}".format(out_file, out_fasta))
    logging.info("aggregated records: {}".format(len(merged)))


def main():
    parser = parse_args()
    args = parser.parse_args()
    if args.datadir == args.outdir:
        logging.info("please specify a different output directory")
        sys.exit(1)

    aggregate_data(datadir=args.datadir, prefix=args.prefix, outdir=args.outdir)


if __name__ == '__main__':
    main()
