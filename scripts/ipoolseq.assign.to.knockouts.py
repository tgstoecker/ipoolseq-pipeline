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

import pysam
import re
import sys
import argparse
import itertools
import copy
import pysamutils
import pprint
from BCBio import GFF
from collections import namedtuple
from collections import defaultdict

# Command-line arguments
parser = argparse.ArgumentParser(description='Assign reads to knockouts')
parser.add_argument('--mapping-fuzzyness', action='store', type=int,
                    default=10, help='Assumed mapping fuzziness')
parser.add_argument('--max-fragment-length', action='store', type=int,
                    default=1000, help='Maximum fragment length (singleton reads only)')
parser.add_argument('input_gff', help='Input GFF2 file with KO cassette insertions', nargs=1)
parser.add_argument('input_bam', help='Input BAM file containing mapped reads', nargs=1)
parser.add_argument('output_bam', help='Output BAM file for assigned reads', nargs=1)
args = parser.parse_args()
ifile_knockouts = args.input_gff[0]
ifile = args.input_bam[0]
ofile = args.output_bam[0]
FEATURE_BOUNDARIES_FUZZ = args.mapping_fuzzyness
FEATURE_FLANK_LENGTH = args.max_fragment_length

print('Reading knockout list from %s' % ifile_knockouts, file=sys.stderr)
print('Reading mapped reads from %s' % ifile, file=sys.stderr)
print('Writing assigned reads to %s' % ofile, file=sys.stderr)
print('Allowed mapping fuzzyness is +/- %d bp, assumed max. fragment length %d bp' % (FEATURE_BOUNDARIES_FUZZ, FEATURE_FLANK_LENGTH), file=sys.stderr)

# Read knockouts from input GFF file
Knockout = namedtuple('Knockout', ['name', 'start', 'end', 'strand', 'left', 'right', 'hits'])
knockouts_bychr = defaultdict(list)
with open(ifile_knockouts) as input_knockouts:
  for r in GFF.parse(input_knockouts):
    knockouts_chr = knockouts_bychr[r.id]
    for f in r.features:
      knockouts_chr.append(Knockout(name=f.qualifiers['Name'][0],
                                  start=int(f.location.start),
                                  end=int(f.location.end),
                                  strand=int(f.strand),
                                  left='5p' if f.strand == +1 else '3p',
                                  right='3p' if f.strand == +1 else '5p',
                                  hits=[0, 0]))

# Sort knockouts within chromosomes by their start position, and
# Build index which allows looking up knockouts by their ID.
knockouts_byid = dict()
for chr, knockouts_chr in knockouts_bychr.items():
  knockouts_chr.sort(key=lambda x: x.start)
  p = None
  for f in knockouts_chr:
    if (p is not None) and (p.end >= f.start):
      raise ValueError('knockouts must not overlap', p.name, f.name)
    p = f
    knockouts_byid[id(f)] = f

# Open input file
input = pysam.AlignmentFile(ifile, mode="rb")
if input.header['HD']['SO'] != 'coordinate':
  raise RuntimeError("input file must be sorted by coordinate")

# Open output file
header = copy.deepcopy(dict(input.header))
output = pysamutils.reorder_buffer(pysam.AlignmentFile(ofile, mode="wb", header=header))

# Information about previously seen mate
Mate = namedtuple('Mate', ['start', 'end'])
mates = {}

fragments = 0
discordant = 0
unmatched = 0
ambiguous = 0
for read in input.fetch(until_eof=True):
  if read.is_unmapped:
    # Skip unmapped reads
    continue
  elif not read.is_proper_pair and not read.mate_is_unmapped:
    # Skip discordant pairs where both mates are mapped
    if read.is_read1:
      discordant = discordant + 1
    continue
  elif read.is_proper_pair and read.query_name not in mates:
    # Wait for mate before processing read
    mates[read.query_name] = (read, output.placeholder())
    continue

  # Switch to correct knockout list if chromosome changed
  r_chr = input.getrname(read.rname)
  if (chr != r_chr):
    knockouts_chr = knockouts_bychr[r_chr]
    chr = r_chr

  # Initialize
  fragments = fragments + 1
  hits = set()

  # Fetch mate (must have been seen earlier)
  start, end = read.reference_start, read.reference_end
  mate, mate_placeholder = None, None
  if read.is_proper_pair:
    mate, mate_placeholder = mates[read.query_name]
    del mates[read.query_name]
    start, end = min(mate.reference_start, start), max(mate.reference_end, end)

  # Find knockouts whose left flank coincides with the fragment's end,
  # or whose right flank coincides with the fragment's start. Knockout boundaries
  # are assumed to be only be determined up to FEATURE_BOUNDARIES_FUZZ bases.
  for f in knockouts_chr:
    if end >= f.start - FEATURE_BOUNDARIES_FUZZ and end < f.start + FEATURE_BOUNDARIES_FUZZ:
      hits.add((id(f), f.left))
    if start >= f.end - FEATURE_BOUNDARIES_FUZZ and start <= f.end + FEATURE_BOUNDARIES_FUZZ:
      hits.add((id(f), f.right))

  # For single-ended reads, additionally consider knockouts a match if the read points
  # toward the knockout and starts no further than FEATURE_FLANK_LENGTH away. Knockout
  # boundaries are again assumed to only be determined up to FEATURE_BOUNDARIES_FUZZ. 
  if not read.is_proper_pair:
    start = read.reference_start
    end = read.reference_end
    if not read.is_reverse:
      for f in knockouts_chr:
        if start >= f.start - FEATURE_BOUNDARIES_FUZZ - FEATURE_FLANK_LENGTH and end <= f.start + FEATURE_BOUNDARIES_FUZZ:
          hits.add((id(f), f.left))
    else:
      for f in knockouts_chr:
        if start >= f.end - FEATURE_BOUNDARIES_FUZZ and end <= f.end + FEATURE_BOUNDARIES_FUZZ + FEATURE_FLANK_LENGTH:
          hits.add((id(f), f.right))
        
  # Determine whether the fragments matches a knockout uniquely
  knockout = None
  if len(hits) == 0:
    unmatched = unmatched + 1
    knockout = 'unmatched'
  elif len(hits) > 1:
    ambiguous = ambiguous + 1
    knockout = 'ambiguous'
  else:
    f_id, s = hits.pop()
    knockout = '%s:%s' % (knockouts_byid[f_id].name, s)

  # Output current read's mate (we process fragments after seeing the second read)
  if mate is not None:
    mate.set_tag(tag='XT', value_type='Z', value=knockout, replace=True)
    output.fillin(mate, placeholder=mate_placeholder)

  # And output this read.
  read.set_tag(tag='XT', value_type='Z', value=knockout, replace=True)
  output.write(read)

output.close()    
print('skipped %d discordant fragments' % discordant, file=sys.stderr)
print('skipped %d fragments matching none of the knockouts' % unmatched, file=sys.stderr)
print('skipped %d fragments matching multiple knockouts' % ambiguous, file=sys.stderr)
