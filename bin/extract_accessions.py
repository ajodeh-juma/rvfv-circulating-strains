#!/usr/bin/env python3
"""
extract sequences from alignment given accessions/sequence headers 
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
    required_group.add_argument('--alignment', required=True, type=str,
                                dest="alignment", metavar="<str>",
                                help="path to the alignment file FASTA format"
                                )
    required_group.add_argument('--metadata', required=True, type=str,
                                dest="metadata", metavar="<str>",
                                help="path to the file containing metadata"
                                )
    required_group.add_argument('--reference', metavar="<str>", dest="reference", default='',
                                help="Reference sequence accession. e.g NC_014398")
    required_group.add_argument('--sort-by-column', required=True, type=str,
                                dest="sort_by_column", metavar="<str>",
                                choices=['accession', 'host', 'country', 'lineage', 'year', 'strain_type'],
                                help="column to sort sequences by (data will be sorted in ascending order) with the "
                                     "first column being the reference sequence",
                                )
    parser.add_argument('--outdir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def extract_sequences(alignment, metadata, reference, sort_by_column, outdir):
    """

    :param alignment:
    :param metadata:
    :param reference
    :param sort_by_column:
    :param outdir
    :return:
    """

    # get country names
    alpha2ctry = dict()
    t = list(pycountry.countries)
    for country in t:
        alpha2ctry[country.alpha_3] = '_'.join(country.name.split(',')[0].split())

    # output
    # outdir = os.path.abspath(os.path.dirname(alignment))

    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    # bs = os.path.splitext(os.path.basename(alignment))[0]
    sample = sort_by_column + '.' + reference + '.sorted.alignment'
    out_fasta = os.path.join(outdir, sample) + '.fasta'
    out_accessions = os.path.join(outdir, sample) + '.all.labels.csv'
    out_labels = os.path.join(outdir, sample) + '.csv'

    # read metadata
    df_meta = pd.read_csv(metadata, sep="\t")
    # df_meta['accession'] = df_meta.taxa.str.split("|").str[0]
    df_meta['year'] = df_meta.date.str.split('-').str[0]
    # df_meta['country'] = df_meta['country'].map(alpha2ctry)

    # read lineages
    # df_lineage = pd.read_csv(lineage)

    # # merge and sort by specified column
    # meta = pd.merge(df_meta, df_lineage, left_on="accession", right_on="Query").rename(columns={'Lineage': 'lineage'})
    meta = df_meta[['accession', 'lineage', 'host', 'country', 'year', 'strain', 'strain_type']]
    # # meta = meta.sort_values(sort_by_column, ascending=True)
    meta = meta.sort_values([sort_by_column, 'lineage', 'year'], ascending=[True, True, True])

    labels = meta[['accession', sort_by_column]].rename(columns={"accession": "name", sort_by_column: "label"})
    # print(meta)

    snipit_labels = meta.copy()
    # snipit_labels["label"] = meta["accession"].str.cat(meta["host"], sep="|").str.cat(meta['lineage'],
    # sep="|").str.cat( meta['year'], sep="|").str.cat(meta['country'], sep="|")
    snipit_labels["label"] = meta["accession"].str.cat(meta["host"], sep="|").str.cat(meta['lineage'], sep="|").str.cat(
        meta[sort_by_column], sep="|").str.cat(meta['year'], sep="|")
    snipit_labels = snipit_labels[['accession', 'label']].rename(columns={"accession": "name"})

    # iterate through the input fasta file
    fasta_dict = dict(fasta_iterator(alignment))

    f_dict = dict()
    for k, v in fasta_dict.items():
        f_dict[k.split("|")[0]] = v

    refs = [reference]
    logging.info("Reference accession: {}".format(''.join(refs)))

    # refs = ['NC_014395', 'NC_014396', 'NC_014397', 'DQ380206']
    reorder_list = list()
    accessions = list(meta.accession)
    # if len(accessions) == 1 and [acc in refs for acc in accessions]:
    #     reorder_list = refs + accessions
    # else:
    for acc in accessions:
        if acc in refs:
            reorder_list.append(acc)
            accessions.remove(acc)
            reorder_list.extend(accessions)
    # print(accessions)
    # print(refs)
    # print(reorder_list)
    logging.info("accessions to extract from alignment {}".format(','.join(reorder_list)))

    with open(out_fasta, 'w') as f_obj:
        for acc in reorder_list:
            if f_dict.get(acc) is None:
                pass
            else:
                f_obj.write(">" + acc + "\n")
                f_obj.write(textwrap.fill(f_dict.get(acc), 80))
                f_obj.write("\n")

    # with open(out_accessions, "w") as f_obj:
    #     f_obj.write("accession" + "\n")
    #     for acc in reorder_list:
    #         f_obj.write(acc + "\n")

    labels.to_csv(out_labels, sep=",", index=False)
    snipit_labels.to_csv(out_accessions, sep=",", index=False)
    return out_fasta, out_accessions, out_labels


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.alignment is not None and isinstance(args.alignment, str):
        args.alignment = args.alignment
        logging.info("input fasta = {:^10}".format(os.path.basename(args.alignment), str))

    if args.metadata is not None and isinstance(args.metadata, str):
        args.metadata = args.metadata
        logging.info("input metadata file = {:^10}".format(os.path.basename(args.metadata), str))

    logging.info("extracting accessions, sorting by {}".format(args.sort_by_column))
    extract_sequences(alignment=args.alignment,
                      metadata=args.metadata,
                      reference=args.reference,
                      sort_by_column=args.sort_by_column,
                      outdir=args.outdir)


if __name__ == '__main__':
    main()
