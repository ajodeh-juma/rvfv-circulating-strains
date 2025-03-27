#!/usr/bin/env python3
"""
convert dna to amino acid
"""

import os
import logging
import argparse
import textwrap
from textwrap import dedent
from Bio import AlignIO

log_level = logging.DEBUG
logging.basicConfig(level=log_level, format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p')


def parse_args():
    # arguments parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        argument_default=argparse.SUPPRESS, prefix_chars='--',
        description=dedent('''convert DNA to protein sequences''')
    )
    # helpstr = """python3 fetch_segments.py [options]"""

    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument('--alignment', required=True, type=str,
                                dest="alignment", metavar="<str>",
                                help="path to the alignment file in fasta format"
                                )
    required_group.add_argument('--outfile', metavar="<str>", dest="outfile", default='',
                                help="path to the output file")
    return parser


def read_alignment(alignment):
    """

    :param alignment:
    :return:
    """

    logging.info("1. Reading in the alignment: {}".format(alignment))
    try:
        aln = AlignIO.read(alignment, "fasta")
    except Exception as error:
        raise AlignmentError("\nERROR: Problem reading in {}: {}".format(alignment, str(error)))
    return aln


def translate_dna(alignment, outfile):
    """
    :param alignment (str) a DNA sequence string
    :return: (str) a protein string from the forward reading frame 1
    """

    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
                  '---': '-',
                  }

    aln = read_alignment(alignment=alignment)

    with open(outfile, 'w') as f_obj:
        for record in aln:
            seq = record.seq.upper()
            prot = []
            for n in range(0, len(seq), 3):
                if seq[n:n + 3] in codontable:
                    residue = codontable[seq[n:n + 3]]
                else:
                    residue = "X"

                prot.append(residue)
            aa = "".join(prot)
            f_obj.write(">" + record.id + "\n")
            f_obj.write(textwrap.fill(str(aa), 80))
            f_obj.write("\n")
    return outfile


def main():
    parser = parse_args()
    args = parser.parse_args()

    if args.alignment is not None and isinstance(args.alignment, str):
        logging.info("input alignment = {:^10}".format(os.path.basename(args.alignment), str))

    translate_dna(alignment=args.alignment, outfile=args.outfile)


if __name__ == '__main__':
    main()
