#! /usr/bin/env python

from __future__ import division
import os
import re
import sys
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

import time
import argparse
from os import path
from Bio import SeqIO
from utils import mkdir



def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

# def splitter(info_file):
#     """
#     split a large file: default split size/chunk size is 1000 lines
#     parameters:
#         - info_file : input file tab delimited format
#     returns:
#         - chunks: list of split files
#     """

#     split_size = 500000
#     file_number = 1
#     file_basename = os.path.splitext(info_file)[0]
#     file_ext = os.path.splitext(info_file)[1]

#     print ("Splitting %s into multiple files with %s lines" % (os.path.join(info_file), str(split_size)))

#     chunks = []
#     keep_headers = True
#     outfile = os.path.join(file_basename)+'_'+str(file_number)+file_ext
#     current_out_writer = csv.writer(open(outfile, 'w'), delimiter=',')
#     current_limit = split_size

#     with open(info_file) as f:
#         reader = csv.reader(f)
#         if keep_headers:
#             headers = next(reader)
#             chunks.append(outfile)
#             current_out_writer.writerow(headers)
#             sys.stdout.write("creating file %s" % (outfile)+'\n')

#         for i, row in enumerate(reader):
#             if i + 1 > current_limit:
#                 file_number += 1
#                 current_limit = split_size * file_number
#                 outfile = os.path.join(file_basename)+'_'+str(file_number)+file_ext
#                 current_out_writer = csv.writer(open(outfile, 'w'), delimiter=',')
#                 if keep_headers:
#                     chunks.append(outfile)
#                     current_out_writer.writerow(headers)
#             current_out_writer.writerow(row)
#         return chunks

# def insert_snp(flanking_fasta, vcf_file):
#     """Function that takes in a flanking file (FASTA format) and a VCF file \n
#     (having SNP chromosome, positions, reference and mutation), modifies the \n
#     flanking FASTA file to include the mutation (alternative)"""
#     vcf_dict = {}
#     fasta_dict = {}
#     the_dict = {}
#     with open(vcf_file) as f_obj:
#         #next(f_obj)
#         for line in f_obj:
#             if 'po' in line:
#                 continue
#             else:
#                 line = line.strip()
#                 vcf_dict[line.split()[0]+' '+(line.split()[1])] = line

#     with open(flanking_fasta) as f_obj:
#         for line in f_obj:
#             if line.startswith('>'):
#                 line = line.strip()
#                 chrm = line.split(':')[0].split('>')[1]
#                 start = line.split(':')[1].split('-')[0]
#                 end = line.split(':')[1].split('-')[1]
#                 key = chrm+' '+start+' '+end
#                 fasta_dict[key] = ''
#             else:
#                 fasta_dict[key] += ''.join(line.strip())
#     out = StringIO()
#     outfile = os.path.splitext(os.path.basename(vcf_file))[0] + '_flanking_sites.txt'
#     outfilename = os.path.join(os.path.dirname(vcf_file), outfile)
#     print ("started writing to {}".format(outfile))
#     with open(outfile, 'w') as f_obj:
#         for k in sorted(fasta_dict.iterkeys()):
#             for key in sorted(vcf_dict.iterkeys()):
#                 if int(k.split()[1]) <= int(key.split()[1]) <= int(k.split()[2]) and (int(key.split()[1]) - int(k.split()[1])) -1 == 21:
#                     fasta = list(fasta_dict[k].lower())
#                     fasta[(int(key.split()[1]) - int(k.split()[1])) -1] = '['+fasta_dict[k][21].upper()+'/'+vcf_dict[key].split()[3].upper()+']'
#                 print ("writing {}\t{}".format(key, ''.join(fasta)))
#                 f_obj.write(key+'\t'+''.join(fasta))
#                 f_obj.write('\n')
#     #fas = out.getvalue()
#     #with open(outfilename, 'w') as fhandle:
#     #    fhandle.write(fas)



parser=argparse.ArgumentParser()
helpstr = """python split_fasta.py [options]"""
parser.add_argument('-f', '--fasta', help="FASTA file")
parser.add_argument('-n', '--num', type=int, help="number of files to generate after splittting")
parser.add_argument('-o', '--outDir', metavar='<DIR>', dest="outdir", help="path to the output directory",
                        default=".")
args=parser.parse_args()

# open input file:
if args.fasta != None:
    fastafile = args.fasta
else:
    sys.stderr.write("Please specify the FASTA file!\n")
    sys.exit(2)

 # create output directory if not exists
outdir = os.path.abspath(args.outdir)
mkdir(outdir)


fastafiles = []
record_iter = SeqIO.parse(open(fastafile),"fasta")
for i, batch in enumerate(batch_iterator(record_iter, args.num)):
    base = os.path.splitext(os.path.basename(fastafile))[0]
    filename = "{}_{}.fasta".format(base, i + 1)
    filename = os.path.join(outdir, filename)
    print(filename)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))
    fastafiles.append(filename)

#r = insert_snp(fastafile, vcffile)

# split = splitter(vcffile)
# for i, j in zip(fastafiles, split):
#     k = insert_snp(i, j)

# # set end time format
# end_time = datetime.now()
# e = end_time.strftime(format)
# tdelta = end_time - start_time

# print("\nCompleted analysis on: "+e+ "\n" )

# # format the time delta object to human readable form
# d = dict(days=tdelta.days)
# d['hrs'], rem = divmod(tdelta.seconds, 3600)
# d['min'], d['sec'] = divmod(rem, 60)

# if d['min'] is 0:
#     fmt = '{sec} sec'
# elif d['hrs'] is 0:
#     fmt = '{min} min {sec} sec'
# elif d['days'] is 0:
#     fmt = '{hrs} hr(s) {min} min {sec} sec'
# else:
#     fmt = '{days} day(s) {hrs} hr(s) {min} min {sec} sec'
# print("[ALL done] Runtime: " +'\t'+fmt.format(**d))
