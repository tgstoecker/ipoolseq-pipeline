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
import pysamutils
import gzip
import pysam
from pysamutils import read_pair_iterator
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

def get_isite_data(read1, read2, bam):
  # Extract chromosome name
  seq = read2.reference_name
  seq_len = bam.get_reference_length(read2.reference_name)

  # Read 2 starts within the inserted transposon and extends outwards into the genomic regions
  # It should thus start (in sequencing direction) with the insertion motif (usually TTAA), and
  # we consider the position of the first base in that motif to be the position of the insertion site.
  if not read2.is_reverse:
    pos = max(read2.reference_start - unmapped_lengths(read2)[0], 0)
    motif = read2.query_sequence[:len(args.motif)]
    dir = "+"
  else:
    pos = min(read2.reference_end + unmapped_lengths(read2)[1], seq_len) - len(args.motif)
    motif = read2.query_sequence[-len(args.motif):]
    dir = "-"

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

# Insertion site sets, indexed by chromosome name and flank
isites_raw = defaultdict(lambda: defaultdict(lambda: set()))

# Find insertion sites
for flank in ['5p', '3p']:
  print("Finding insertion sites for %s flank" % flank, file=sys.stderr)

  ifile = {"5p": args.input_bam_5p[0], "3p": args.input_bam_3p[0]}[flank]
  input = pysam.AlignmentFile(ifile, mode="rb")
  for read1, read2 in read_pair_iterator(input, item=[read_pair_iterator.READ_1, read_pair_iterator.READ_2]):
    # Get insertion site - related information
    seq, pos, dir, motif = get_isite_data(read1, read2, input)

    # Ignore reads where the motif differes from the expected motif by more than a single base
    if (hamming_distance(motif, args.motif) > 1):
      continue

    # Store insertion site
    isites_raw[seq]["%s%s" % (flank, dir)].add(pos)

  # Done scanning reads
  input.close()
  
# Filter insertion sites. Both 5' and 3' fragments must be found for a site,
# and the two must map in opposite directions
isites = {}
for seq in isites_raw.keys():
  isites[seq] = ((isites_raw[seq]["5p+"] & isites_raw[seq]["3p-"]) |
                 (isites_raw[seq]["5p-"] & isites_raw[seq]["3p+"]))
print("Found %d insertion sites with 3p and 5p reads" % sum(len(s) for s in isites.values()), file=sys.stderr)

# Output valid insertion sites as GFF3
if args.output_gff is not None:
  gff = gzip.open(args.output_gff, 'wt', compresslevel=6)
  for seq in isites.keys():
    for pos in isites[seq]:
      gff.write("{seq}\t.\t.\t{start}\t{end}\t{score}\t{strand}\t.\t{attrs}\n".format(
                seq=seq, start=pos+1, end=pos+len(args.motif), score=".", strand=".", attrs=""))
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
    seq, pos, dir, motif = get_isite_data(read1, read2, input)

    # Set XT tag to insertion site
    if pos in isites[seq]:
      xt = "%s|%d:%s" % (seq, pos, flank)
      read1.set_tag(tag="XT", value_type="Z", value=xt)
      read2.set_tag(tag="XT", value_type="Z", value=xt)
    else:
      read1.set_tag(tag="XT", value_type="Z", value="unassigned")
      read2.set_tag(tag="XT", value_type="Z", value="unassigned")

    read_pairs.write(read1, read2)

  # Done with the flank
  output.close()    
