#!/usr/bin/env python3
"""
summarize SNPs per grouping
"""

import os
import re
import csv
import logging
import argparse
import textwrap
import pandas as pd
from textwrap import dedent
from collections import defaultdict
from collections import Counter


from utils import str2bool
from utils import mkdir
from utils import read_alignment

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent('''get SNP records''')
    )
    # helpstr = """python3 fetch_segments.py [options]"""

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--alignment', required=True, type=str,
                                dest="alignment", metavar="<str>",
                                help="path to the alignment file in fasta format"
                                )
    required_group.add_argument('--reference', metavar="<str>", dest="reference", default='',
                                help="Reference sequence accession. e.g NC_014398")
    required_group.add_argument('--metadata', required=True, type=str,
                                dest="metadata", metavar="<str>",
                                help="path to the file containing metadata"
                                )
    required_group.add_argument('--min-freq', type=float,
                                dest="min_freq", metavar="<float>", default=0.75,
                                help="proportion of the minimum number of sequences which contain the SNP"
                                )
    parser.add_argument('--seq-type', type=str, metavar="<str>", dest="seq_type", default='dna',
                        choices=['dna', 'protein'],
                        help="specify the type of sequence either 'dna' or 'protein'"
                        )
    required_group.add_argument('--max-freq', type=float,
                                dest="max_freq", metavar="<float>", default=1.0,
                                help="proportion of the maximum number of sequences which contain the SNP"
                                )
    required_group.add_argument('--snp-type', required=True, type=str,
                                dest="snp_type", metavar="<str>", choices=['singleton', 'multiple', 'conserved', 'all'],
                                default="all",
                                help="snp type, singleton (only a single snp per position), multiple (more than one "
                                     "snp per position) and conserved (snps that occur commonly across all the "
                                     "sequences)",
                                )
    required_group.add_argument('--prefix', required=True, metavar='<str>',
                                dest="prefix",
                                help="prefix of the output file(s)"
                                )
    parser.add_argument('--outdir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
    return parser




def get_reference(alignment, reference_accession):
    """
    get the reference sequence from an alignment given the reference id

    :param alignment:
    :param reference_accession:
    :return:
    """

    aln = read_alignment(alignment=alignment)
    for record in aln:
        if record.id.split('|')[0] == reference_accession:
            return record


def find_snps(reference, sequence, seq_type):
    """
    :param reference:
    :param sequence:
    :return:
    """




    nucleotides = ['A', 'C', 'T', 'G', '-']
    amino_acids = ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-"]

    snps = list()
    index = 0

    if seq_type == 'dna':
        for i in range(len(reference)):
            if reference[i] != 'N':
                index += 1
            col = [reference[i], sequence[i]]
            if len(set(col)) > 1:
                if not col[1].upper() in nucleotides:
                    pass
                else:
                    snp = f"{col[0].upper()}{index}{col[1].upper()}"
                    snps.append(snp)
    else:
        for i in range(len(reference)):
            if reference[i] != 'X':
                index += 1
            col = [reference[i], sequence[i]]
            if len(set(col)) > 1:
                if not col[1].upper() in amino_acids:
                    pass
                else:
                    snp = f"{col[0].upper()}{index}{col[1].upper()}"
                    snps.append(snp)

    return snps

def pcent_done(c, total):
    """

    :param c:
    :param total:
    :return:
    """
    return round((c * 100) / total, 2)



def snp_list_to_snp_string(snp_list):
    """
    convert a snp list into a ; separated string of snps sorted by position
    :param snp_list:
    :return:
    """

    # snp_string = ";".join(sorted(snp_list, key=lambda x: int(x[:-2])))
    snp_string = ";".join(sorted(snp_list, key=lambda x: int(x[1:-1])))
    return snp_string

def append_snps_annotations(alignment, reference, seq_type):
    """
    add list of snps relative to ref as an annotation to the seq record
    :param alignment:
    :param reference:
    :param seq_type:
    :return:
    """

    d = defaultdict()

    reference = get_reference(alignment=alignment, reference_accession=reference)
    alignment = read_alignment(alignment=alignment)

    c = 0
    total = len(alignment)
    for record in alignment:
        c += 1
        if c % 500 == 0:
            print(pcent_done(c, total), '%')
        snps = find_snps(reference.seq, record.seq, seq_type)
        record.annotations["snps"] = snps
        snp_string = snp_list_to_snp_string(snps)
        record.annotations["snp_string"] = snp_string

        d[record.id] = record.annotations['snp_string']
    logging.info("{} records annotated".format(total))
    return d




def write_snp_positions(snps_records_dict, metadata, min_freq, max_freq, snp_type, prefix, outdir):
    """

    :param snps_records_dict:
    :param metadata:
    :param min_freq:
    :param max_freq:
    :param snp_type:
    :param prefix:
    :param outdir:
    :return:
    """

    outdir = os.path.abspath(outdir)
    mkdir(outdir)

    out1 = os.path.join(outdir, prefix) + '.mutations.per.strain.' + snp_type + '.csv'
    out2 = os.path.join(outdir, prefix) + '.' + snp_type + '.txt'


    df = pd.DataFrame(snps_records_dict.items(), columns=["taxa", "snps"])
    df['positions'] = df['snps'].apply(lambda x: re.findall(r"[-+]?\d*\.\d+|\d+", str(x))).map(
        lambda x: [int(i) if '.' not in i else float(i) for i in x]).tolist()
    df['SNPs'] = df['snps'].apply(lambda x: re.findall(r"([A-Z]+)", str(x)))
    df['SNPs'] = df['SNPs'].apply(lambda x: [''.join(t) for t in list(zip(x[::2], x[1::2]))])
    df['accession'] = df.taxa.str.split("|").str[0]
    df.snps = df.snps.str.split(';')

    # read metadata
    df1 = pd.read_csv(metadata, sep="\t")
    df2 = pd.merge(df, df1, on='accession').rename(columns={'taxa_y': 'taxa'})


    # columns of interest
    columns = ['accession', 'taxa', 'lineage', 'strain', 'strain_type', 'host', 'year', 'snps', 'positions']
    df3 = df2[columns]


    with open(out1, 'w') as f_obj:
        writer = csv.writer(f_obj)
        writer.writerow(columns)
        for index, row in df3.iterrows():
            for i in range(len(row['positions'])):
                l = [row['accession'], row['taxa'], row['lineage'], row['strain'], row['strain_type'], row['host'], row['year'],  row['snps'][i], row['positions'][i]]
                writer.writerow(l)


    #
    to_sum = len(df)
    min_seqs = int(min_freq * to_sum)
    max_seqs = int(max_freq * to_sum)

    logging.info("Total number of sequence: {}".format(to_sum))
    logging.info("Minimum number of sequences with {} SNPs: {}".format(snp_type, min_seqs))
    logging.info("Maximum number of sequences with {} SNPs: {}".format(snp_type, max_seqs))


    snp_positions = list()
    snps_counts_per_position = dict()
    for index, row in df.iterrows():
        positions = row['positions']
        snp_positions += [pos for pos in positions]

        for pos in positions:
            if pos not in snps_counts_per_position:
                snps_counts_per_position[pos] = [pos]
            else:
                snps_counts_per_position[pos].append(pos)

    result = list()
    for pos, counts in snps_counts_per_position.items():
        if snp_type == 'singleton':
            if len(counts) > min_seqs:
                result.append(pos)
        elif snp_type == 'conserved':
            if len(counts) >= max_seqs:
                result.append(pos)
        elif snp_type == 'multiple':
            if len(counts) <= min_seqs:
                result.append(pos)
        else:
            result.append(pos)

    with open(out2, 'w') as f_obj:
        f_obj.write(' '.join(str(item) for item in sorted(result)))
    return out1, out2




def main():
    parser = parse_args()
    args = parser.parse_args()

    if args.metadata is not None and isinstance(args.metadata, str):
        logging.info("input metadata = {:^10}".format(os.path.basename(args.metadata), str))

    outdir = os.path.abspath(args.outdir)
    mkdir(outdir)


    d = append_snps_annotations(alignment=args.alignment,
                                reference=args.reference,
                                seq_type=args.seq_type)

    write_snp_positions(snps_records_dict=d,
                        metadata=args.metadata,
                        min_freq=args.min_freq,
                        max_freq=args.max_freq,
                        snp_type=args.snp_type,
                        prefix=args.prefix,
                        outdir=args.outdir)





if __name__ == '__main__':
    main()
