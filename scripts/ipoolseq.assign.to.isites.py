#!/usr/bin/env python
# ipoolseq.assign.to.knockouts.py, Copyright 2016, 2017, 2019 Florian G. Pflug
# 
# This file is part of the iPool-Seq Analysis Pipeline
#
# The iPool-Seq Analysis Pipeline is free software: you can redistribute it
# and/or modify it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# The iPool-Seq Analysis Pipeline is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with the iPool-Seq Analysis Pipeline.  If not, see
# <http://www.gnu.org/licenses/

import re
import sys
import argparse
import copy
import gzip
import urllib.parse
import pysam
import pysamutils
from pysamutils import read_pair_iterator
from itertools import chain, groupby
from collections import defaultdict

# Returns the number of positions in which strings s1 and s2 differ
def hamming_distance(s1, s2):
  """Computed the (unscaled) hamming distance between two lists"""
  if len(s1) != len(s2):
    raise RuntimeError(
        "hamming_distance required sequences with the same length, but %d != %d characters" % (len(s1), len(s2)))
  return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# CIGAR operators that consume query but not reference
CIGAR_QUERY_ONLY = {pysam.CINS: True, pysam.CSOFT_CLIP: True}

# Returns a tuple (prefix, postfix) of the number of unmapped bases and the
# beginning respectively end of the read (in reference direction)
def unmapped_lengths(read):
  if read.is_unmapped:
    return (None, None)
  cigar = read.cigartuples
  pre, post = 0, 0
  for i, (op, len) in enumerate(cigar):
    if op not in CIGAR_QUERY_ONLY:
      break
    pre = pre + len
  for i, (op, len) in enumerate(reversed(cigar)):
    if op not in CIGAR_QUERY_ONLY:
      break
    post = post + len
  return (pre, post) 

def get_isite_data(read1, read2, bam, flank):
  # Extract chromosome name
  seq = read2.reference_name
  seq_len = bam.get_reference_length(read2.reference_name)

  # Read 2 starts within the inserted transposon and extends outwards into the genomic regions
  # It should thus start (in sequencing direction) with the insertion motif (usually TTAA), and
  # we consider the position of the first base in that motif to be the position of the insertion site.
  if not read2.is_reverse:
    pos = max(read2.reference_start - unmapped_lengths(read2)[0], 0)
    motif = read2.query_sequence[:len(args.motif)]
    dir = {"5p": "-", "3p": "+"}[flank]
  else:
    pos = min(read2.reference_end + unmapped_lengths(read2)[1], seq_len) - len(args.motif)
    motif = read2.query_sequence[-len(args.motif):]
    dir = {"5p": "+", "3p": "-"}[flank]

  return (seq, pos, dir, motif)

# Command-line arguments
parser = argparse.ArgumentParser(description='Assign reads to knockouts')
parser.add_argument('--motif', metavar="MOTIF", default="TTAA",
                    help='Expected motif at insertion sites')
parser.add_argument('--output-gff', metavar="GFF3",
                    help='Output GFF3 file with insertion sites')
parser.add_argument('input_bam_5p', help='Input BAM file containing mapped reads for the 5p fragments', nargs=1)
parser.add_argument('input_bam_3p', help='Input BAM file containing mapped reads for the 3p fragments', nargs=1)
parser.add_argument('output_bam_5p', help='Output BAM file for assigned reads of the 5p fragments', nargs=1)
parser.add_argument('output_bam_3p', help='Output BAM file for assigned reads of the 3p fragments', nargs=1)
args = parser.parse_args()

print('Expected motif %s' % args.motif, file=sys.stderr)
print('Reading mapped reads from %s and %s' % (args.input_bam_5p[0], args.input_bam_3p[0]), file=sys.stderr)
print('Writing assigned reads to %s and %s' % (args.output_bam_5p[0], args.output_bam_5p[0]), file=sys.stderr)
if args.output_gff is not None:
  print('Writing insertion sites in GFF3 format to %s' % args.output_gff, file=sys.stderr)

# Insertion site sets
#
# An insertion at a particular site can have two possible orientations, "+" and "-".
# For a "+"-strand insertion, reads originating at the 5' flank of the inserted cassette
# map to the negative strand of the genome (i.e. map in reverse direction) and reads originating
# at the 3' flank map to the positive strand (i.e. in forward direction). For a "-"-strand insertion,
# its exactly the opposite (see also get_isite_data).
#
# Note that "5' flank" and "3' flank" refer to an arbitrarily chosen "reference strand" of the double-
# stranded insertion, which is identified through its sequence (see ipoolseq.trim.py), and note that this
# assumes that the inserted cassette is NOT palindromic.
#
# Tentative insertion sites found in reads are stored in isites_tentative and are indexed by
# a three-level key comprising the chromosome, the strand, and the flank from which the reads
# originated.
#
# Confirmed insertion sites are those for which both 5'- and 3'-originating reads where found,
# these are stored in isites, indexed by chromosone and strand.
isites_tentative = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: set())))

# Find tentative insertion sites
for flank in ['5p', '3p']:
  print("Finding insertion sites for %s flank" % flank, file=sys.stderr)

  ifile = {"5p": args.input_bam_5p[0], "3p": args.input_bam_3p[0]}[flank]
  input = pysam.AlignmentFile(ifile, mode="rb")
  for read1, read2 in read_pair_iterator(input, item=[read_pair_iterator.READ_1, read_pair_iterator.READ_2]):
    # Get insertion site - related information
    seq, pos, dir, motif = get_isite_data(read1, read2, input, flank)

    # Ignore reads where the motif differes from the expected motif by more than a single base
    if (hamming_distance(motif, args.motif) > 1):
      continue

    # Store insertion site
    isites_tentative[seq][dir][flank].add(pos)

  # Done scanning reads
  input.close()
  
# Filter insertion sites by looking for sites which a 5'- and a 3'-signal.
isites = defaultdict(lambda: defaultdict(lambda: set()))
for seq, strands in isites_tentative.items():
  for s in strands.keys():
    isites[seq][s] = isites_tentative[seq][s]["3p"] & isites_tentative[seq][s]["5p"]
print("Found %d +strand insertion sites confirmed by both 3p and 5p reads" % sum(len(strands["+"]) for strands in isites.values()), file=sys.stderr)
print("Found %d -strand insertion sites confirmed by both 3p and 5p reads" % sum(len(strands["-"]) for strands in isites.values()), file=sys.stderr)
print("Found %d insertion sites in total confirmed by both 3p and 5p reads" % sum(len(strands["+"] | strands["-"]) for strands in isites.values()), file=sys.stderr)

# Output valid insertion sites as GFF3
if args.output_gff is not None:
  gff = gzip.open(args.output_gff, 'wt', compresslevel=6)
  for seq, strands in isites.items():
    for pos, group in groupby(sorted(chain(strands["+"], strands["-"]))):
      for s in strands.keys():
        if pos in strands[s]:
          name = "%s:%d(%s)" % (seq, pos, s)
          attrs = { "ID": name, "Name": name }
          gff.write("{seq}\t.\t{type}\t{start}\t{end}\t{score}\t{strand}\t.\t{attrs}\n".format(
                    seq=seq, type="transposable_element_insertion_site",
                    start=pos+1, end=pos+len(args.motif), score=".", strand=s,
                    attrs=";".join("%s=%s" % (urllib.parse.quote(k), urllib.parse.quote(v))
                                   for k, v in attrs.items())))
  gff.close()

# Assign reads
for flank in ['5p', '3p']:
  print("Assigning reads for %s flank" % flank, file=sys.stderr)

  ifile = {"5p": args.input_bam_5p[0], "3p": args.input_bam_3p[0]}[flank]
  ofile = {"5p": args.output_bam_5p[0], "3p": args.output_bam_3p[0]}[flank]

  # Assign reads to valid insertion sites
  input = pysam.AlignmentFile(ifile, mode="rb")
  header = copy.deepcopy(dict(input.header))
  output = pysamutils.reorder_buffer(pysam.AlignmentFile(ofile, mode="wb", header=header))
  read_pairs = read_pair_iterator(input, item=[read_pair_iterator.READ_1, read_pair_iterator.READ_2], output=output)
  for read1, read2 in read_pairs:
    # Get insertion site - related information
    seq, pos, dir, motif = get_isite_data(read1, read2, input, flank)

    # Set XT tag to insertion site
    if seq in isites and pos in isites[seq][dir]:
      xt = "%s:%d(%s)" % (seq, pos, dir)
      read1.set_tag(tag="XT", value_type="Z", value=xt)
      read2.set_tag(tag="XT", value_type="Z", value=xt)
    else:
      read1.set_tag(tag="XT", value_type="Z", value="unassigned")
      read2.set_tag(tag="XT", value_type="Z", value="unassigned")

    read_pairs.write(read1, read2)

  # Done with the flank
  output.close()    
