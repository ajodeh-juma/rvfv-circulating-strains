#!/usr/bin/env python3
"""
summarize SNPs per grouping
"""

import os
import logging
import argparse
from textwrap import dedent
import pandas as pd
from collections import defaultdict

from utils import mkdir
from utils import find_executable
from utils import run_shell_command
from utils import str2bool

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent('''Compute TsTv stats and parse the file''')
    )
    # helpstr = """python3 fetch_segments.py [options]"""

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--vcf', required=True, type=str,
                                dest="vcf", metavar="<str>",
                                help="path to the vcf file having variants"
                                )
    required_group.add_argument('--metadata', required=True, type=str,
                                dest="metadata", metavar="<str>",
                                help="path to the file containing metadata"
                                )
    required_group.add_argument('--reference', metavar="<str>", dest="reference", default='',
                                help="Reference sequence accession. e.g NC_014398")
    parser.add_argument('--recomb', metavar="<str>", dest="recomb", required=False,
                        help="path to a text file having one or more accessions recombinant sequences with a header named 'recombinants'"
                        )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--group-per-lineage', type=str2bool, metavar="<bool>",
                        default=True, const=True, nargs='?',
                        dest='group_per_lineage',
                        help="specify if you want group the stats output per lineage"
                        )
    parser.add_argument('--outdir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser


def run_snpsift_tstv(vcf, prefix, outdir):
    """

    :param vcf:
    :param prefix:
    :param outdir:
    :return:
    """

    # create output directory
    outdir = os.path.abspath(outdir)
    mkdir(outdir)
    outfile = os.path.join(outdir, prefix + '.snpSift.tstv.txt')

    # locate the executable
    exe = find_executable(['SnpSift'])

    # if os.path.exists(outfile):
    #     logging.critical("SnpSift TsTv file {} exists!".format(os.path.basename(outfile)))
    # else:
    call = ["{} tstv {} > {}".format(exe, vcf, outfile)]
    cmd = " ".join(call)
    run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return outfile


def parse_snpsift_tstv(snpsift_tstv, metadata, reference, recomb, prefix, group_per_lineage, outdir):
    """

    :param snpsift_tstv:
    :param metadata:
    :param reference:
    :param recomb:
    :param prefix:
    :param group_per_lineage:
    :param outdir:
    :return:
    """

    df = pd.read_csv(recomb, sep="\t")
    print(df)
    logging.info("Recombinants: {}".format(list(df.recombinants)))

    refs = [reference]
    logging.info("Reference accession: {}".format(''.join(refs)))

    # create output directory
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    outfile = os.path.join(outdir, prefix + '.all.samples.snps.csv')

    if group_per_lineage:

        df1 = pd.read_csv(metadata, sep="\t")
        df2 = pd.read_csv(snpsift_tstv, sep=",")
        # df2 = df2.drop(4)

        df2 = df2.reset_index().T

        df2 = df2[[0, 1, 2, 3]].drop(['level_0', 'TS/TV stats:']).rename(
            columns={0: 'accession', 1: 'Ti', 2: 'Tv', 3: 'Ti/Tv'}).reset_index().drop(['index'], axis=1)

        df3 = pd.merge(df1, df2, on='accession')
        df3 = df3.round(1)

        # remove reference accession and recombinants
        df3 = df3[~df3['accession'].isin(refs)]

        if len(list(df.recombinants)) > 0:
            df3 = df3[~df3['accession'].isin(list(df.recombinants))]

        df3.to_csv(outfile, index=False, header=True)

    return outfile


def run_bcftools_stats(vcf, prefix, outdir):
    """

    :param vcf:
    :param prefix:
    :param outdir:
    :return:
    """

    # create output directory
    outdir = os.path.abspath(outdir)
    mkdir(outdir)
    outfile = os.path.join(outdir, prefix + '.bcftools.stats.tstv.txt')

    # locate the executable
    exe = find_executable(['bcftools'])

    # if os.path.exists(outfile):
    #     logging.critical("bcftools TsTv file {} exists!".format(os.path.basename(outfile)))
    # else:
    call = ["{} stats -s - {} > {}".format(exe, vcf, outfile)]
    cmd = " ".join(call)
    run_shell_command(cmd=cmd, raise_errors=False, extra_env=None)
    return outfile


def parse_bcftools_stats(bcftools_stats, prefix, outdir):
    """

    :param bcftools_stats:
    :param prefix:
    :param outdir:
    :return:
    """

    sn_dict = defaultdict(int)
    tstv_dict = defaultdict()
    st_dict = defaultdict(int)
    with open(bcftools_stats) as f_input:
        for line in f_input:
            if line.startswith("#"):
                continue
            else:
                line = line.strip()
                if line.startswith("SN"):
                    # samples, records, no-ALTs, SNPs, MNPs, INDELs, others, Multiallelic
                    key = line.split('\t')[2].rsplit('number of ')[1].rsplit(':')[0]
                    value = int(line.split()[-1])
                    if key not in sn_dict:
                        sn_dict[key] = value

                if line.startswith("TSTV"):
                    # ts, tv, tstv, ts(1st ALT), tv(1st ALT), ts/tv(1st ALT)
                    tstv_dict['ts'] = int(line.split("\t")[2])
                    tstv_dict['tv'] = int(line.split("\t")[3])
                    tstv_dict['tstv'] = float(line.split("\t")[4])
                    tstv_dict['tsALT'] = int(line.split("\t")[5])
                    tstv_dict['tvALT'] = int(line.split("\t")[6])
                    tstv_dict['tstvALT'] = float(line.split("\t")[7])

                if line.startswith("ST"):
                    key = line.split("\t")[2]
                    value = int(line.split("\t")[-1])
                    # A>C, A>G, A>T, C>A, C>G, C>T, G>A, G>C, G>T, T>A, T>C, T>G
                    if key not in st_dict:
                        st_dict[key] = value

    df1 = pd.DataFrame(sn_dict.items()).T
    df2 = pd.DataFrame(tstv_dict.items()).T
    df3 = pd.DataFrame(st_dict.items()).T
    df = pd.concat([df1, df2, df3], axis=1)

    # create output directory
    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    outfile = os.path.join(outdir, prefix + '.parsed.bcftools.stats.csv')
    df.to_csv(outfile, index=False, header=False)
    return outfile


def main():
    parser = parse_args()
    args = parser.parse_args()

    # check the arguments
    if args.vcf is not None and isinstance(args.vcf, str):
        logging.info("input vcf = {:^10}".format(os.path.basename(args.vcf), str))

    if args.metadata is not None and isinstance(args.metadata, str):
        logging.info("input metadata = {:^10}".format(os.path.basename(args.metadata), str))

    bcftools_stats = run_bcftools_stats(vcf=args.vcf, prefix=args.prefix, outdir=args.outdir)
    parse_bcftools_stats(bcftools_stats=bcftools_stats,
                         prefix=args.prefix,
                         outdir=args.outdir
                         )

    snpsift_tstv = run_snpsift_tstv(vcf=args.vcf, prefix=args.prefix, outdir=args.outdir)
    parse_snpsift_tstv(snpsift_tstv=snpsift_tstv,
                       metadata=args.metadata,
                       reference=args.reference,
                       recomb=args.recomb,
                       prefix=args.prefix,
                       group_per_lineage=args.group_per_lineage,
                       outdir=args.outdir)


if __name__ == '__main__':
    main()
