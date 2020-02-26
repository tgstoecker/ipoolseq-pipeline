# pysamutils.py, Copyright 2016, 2017, 2019 Florian G. Pflug
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

import sys
import pysam
from collections import deque
from collections import defaultdict
from namedlist import namedlist


class genome_cache:
  """
  Cached access of genomes in FASTA format

  When a particular region is requested, a bigger chunk (1 Mb by default)
  starting at the requested region is fetched and cached. Later region
  requests are fulfilled from the cache if possible.
  """

  def __init__(self, filename, chunksize=1024*1024):
    self.chunksize = chunksize
    self.genome = pysam.FastaFile(filename)
    self.curr_seq = None
    self.curr_start = None
    self.curr = None

  def length(self, seq):
    return self.genome.get_reference_length(seq)

  def get(self, seq, start, end):
    if self.curr_seq != seq or self.curr_start > start or self.curr_start + len(self.curr) < end:
      # Requested range does not overlap currently cached region
      if (end > self.genome.get_reference_length(seq)):
        raise IndexError("Requested region [%d, %d) exceeds length of %s (%d)" %
                         (start, end, seq, self.genome.get_reference_length(seq)))
      reg = self.genome.fetch(
          seq, start, start + max(start - end, self.chunksize))
      self.curr_seq, self.curr_start, self.curr = seq, start, reg
    # Return requested region
    return self.curr[(start - self.curr_start):(end - self.curr_start)]


class read_pair_iterator:
  """
  Iterate over the read in a BAM file as mate pairs

  Reads are generally reported in the order they appear in the input file,
  which must be sorted by genomic coordinate (and which unmapped reads last).
  Reads that a part of a proper pair (i.e. flagged as such in the BAM file) are
  reported as pairs, with the leftmost (in genomic coordinates) mate of a pair
  deciding the reporting order.

  Mate pairs where only one mate is mapped are reported as singleton reads if
  skip_unpaired is set to false, or skipped if set to true.

  Mate pairs where both mates are mapped, but in a discordant fasion (i.e. not
  flagged as proper in the BAM file, meaning their distance was too large, they
  ended up on different chromosomes, or had the wrong orientations) are skipped.

  For each pair a tuple with a configurable layout (via the constructor"s item
  parameter) is returned. The tuple can consist of arbitrary combinations of

    * `READ_1`, meaning the 1st read in the pair (according to the sequencer)
    * `READ_2`, meaning the 2nd read in the pair (according to the sequencer)
    * `READ_FWD`, meaning the forward-mapping read in the pair
    * `READ_REV`, meaning the reverse-mapping read in the pair
    * `READ_LEFT`, meaning the leftmost read (in genomic coordinates)
    * `READ_RIGHT`, meaning the rightmost read (in genomic coordinates)
  """

  # Queue entries
  read_pair = namedlist("read_pair", ['left', 'right', 'left_placeholder', 'right_placeholder', 'complete'])

  # Possible entries in the <item> parameter
  READ_1 = 0
  READ_2 = 1
  READ_FWD = 2
  READ_REV = 3
  READ_LEFT = 4
  READ_RIGHT = 5

  # Returns item entry for a given entry type
  item_entry = {
    READ_1: lambda left, right: 0 if left.is_read1 else 1,
    READ_2: lambda left, right: 0 if left.is_read2 else 1,
    READ_FWD: lambda left, right: 0 if not left.is_reverse else 1,
    READ_REV: lambda left, right: 0 if left.is_reverse else 1,
    READ_LEFT: lambda left, right: 0,
    READ_RIGHT: lambda left, right: 1
  }

  def __init__(self, read_iterator, item=(READ_FWD, READ_REV), output=None, skip_unpaired=True):
    # Settings
    self.skip_unpaired = skip_unpaired
    self.item = item
    # Input iterator
    self.iterator = read_iterator
    # Queue and read name index into queue
    self.sequence = None
    self.queue = deque()
    self.left_names = dict()
    # Current item
    self.current_pair = None
    self.current_item_read_indices = None
    self.current_item = None
    self.current_item_written = False
    # Output queue (must be an instance of reorder_buffer)
    self.output = output
    # Statistics
    self.skipped = defaultdict(int)
    self.total = 0

  def __iter__(self):
    return self

  def __next__(self):
    # Repeat until a valid pair is found
    while True:
      # I. If the previous pair wasn't written to the output, remove the reads placeholders
      if self.output is not None and self.current_pair is not None and not self.current_item_written:
        self.output.dontfill(self.current_pair.left_placeholder)
        self.output.dontfill(self.current_pair.right_placeholder)
        self.current_pair.left_placeholder = None
        self.current_pair.right_placeholder = None
      self.current_item_written = False

      # II. Consume input until a complete pair is available
      #
      # Scan ahead until the front of the queue contains a complete pair.
      while not self.queue or not self.queue[0].complete:
        # Queue is empty or first entry is not yet complete.

        # 1. Fetch read
        #
        try:
          r = self.iterator.__next__()
        except StopIteration:
          # Sequences changes to 'None'
          self.on_sequence_will_change()
          raise
        # Check if sequence name changes, run checks if it does
        if self.sequence is not None and self.sequence != r.reference_name:
          self.on_sequence_will_change()
        self.sequence = r.reference_name

        # 2. Enqueue or skip
        #
        # (A) If both mates are unmapped, the pair is ignored. (B) If exactly
        # one mate is mapped, the pair is accepted, and the other mate
        # set to "None". If skip_unpaired was set, the pair is dropped instead.
        # (C) If both mates are mapped but discordantely, the pair is ignored.
        # (D) If both are mapped concordantly, the pair is (naturally) accepted.
        if r.is_unmapped and r.mate_is_unmapped:
          # A. Ignore.
          if r.is_read1:
            self.total += 1
            self.skipped['unmapped'] += 1
        elif r.is_unmapped and not r.mate_is_unmapped:
          # B, unmapped mate. Ignore
          pass
        elif r.mate_is_unmapped and not r.is_unmapped:
          self.total += 1
          if not self.skip_unpaired:
            # B, mapped mate. Create pair entry with missing mate, enqueue.
            r_ph = self.output.placeholder() if self.output is not None else None
            p = read_pair_iterator.read_pair(r, None, r_ph, None, True)
            self.queue.append(p)
          else:
            # B, mapped mate. Ignore because skip_unpaired is set
            self.skipped['unpaired'] += 1
        elif not r.is_proper_pair:
          # C. Ignore.
          if r.is_read1:
            self.total += 1
            self.skipped['discordant'] += 1
        elif r.query_name in self.left_names:
          # D, and mate was already seen, i.e. pair already has an entry.
          p = self.left_names[r.query_name]
          del self.left_names[r.query_name]
          p.right = r
          p.right_placeholder = self.output.placeholder() if self.output is not None else None
          p.complete = True
        else:
          # D, and pair has no entry yet. Create one and flag as incomplete.
          r_ph = self.output.placeholder() if self.output is not None else None
          p = read_pair_iterator.read_pair(r, None, r_ph, None, False)
          self.queue.append(p)
          self.left_names[r.query_name] = p
          self.total += 1

      # III. Dequeue pair and return
      #
      # Dequeue and return frontmost pair, which is now known to be complete
      # We assume that exactly one of the reads maps to the forward strand,
      # and return that one as the first read. If the forward read lies fully
      # to the right of the reverse read (i.e., if they look like this:
      # <--| |-->), the pair is skipped.
      p = self.queue.popleft()
      self.current_pair = p
      # Validate current pair
      if (p.right is not None) and (p.left.is_reverse == p.right.is_reverse):
        raise RuntimeError(
            "concordant pairs expected to have orientation FR or RF")
      if ((p.left is not None) and (p.right is not None) and
          (p.left.reference_end <= p.right.reference_start) and p.left.is_reverse):
        self.skipped['outwards'] += 1
        continue
      # Construct tuple to return, according to the layout desired by the user (self.item)
      # We store the indices (0 for left, 1 for right) of all reads in the returned item
      # in self.current_item_read_indices, which allows write() to figure out which reads
      # was the "left" and which read was the "right" read.
      self.current_item_read_indices = [ self.item_entry[e](p.left, p.right) for e in self.item ]
      self.current_item = [ p.left if i == 0 else p.right for i in self.current_item_read_indices ]
      return self.current_item

  def write(self, read_a, read_b):
    if self.output is None:
      return
    if self.current_item_written:
      raise RuntimeError("write() may be invoked at most once per input pair")
    # Lookup output placeholders for the reads and fill in reads instead
    try:
      for i in range(0, 2):
        if self.current_item_read_indices[i] == 0:
          ph = self.current_pair.left_placeholder
          self.current_pair.left_placeholder = None
        elif self.current_item_read_indices[i] == 1:
          ph = self.current_pair.right_placeholder
          self.current_pair.right_placeholder = None
        else:
          raise RuntimeError("invalid index in current_item_read_indices")
        self.output.fillin(read_a if i == 0 else read_b, ph)
    finally:
      self.current_item_written = True

  # Iterator reached the next sequence.
  # No proper pairs should extend over multiple sequences, so warn if that is the case
  def on_sequence_will_change(self):
    incomplete = sum(not e.complete for e in self.queue)
    if self.queue:
      print("WARNING: %d incomplete concordant mate pairs left af the end of %s" %
            (incomplete, self.sequence),
            file=sys.stderr)
    if incomplete != len(self.left_names):
      raise RuntimeError("inconsistent queue state, %d incomplete pairs but %d left_names entries" %
                         (incomplete, len(self.left_names)))
    if (self.output is not None) and (self.output.nr_postponed() > 0):
      raise RuntimeError("inconsistent output state, %d postponed reads at the end of %s" %
                         (self.output.nr_postponed(), self.sequence))
    self.queue.clear()
    self.left_names.clear()


class merge_reads_iterator:
  """
  Iterate over multiple coordinate-sorted BAM files simulatenously

  The reads from the individual files are merged-joined.

  *NOTE:* The code does not check whether read names are duplicated between
  files. For Illumina reads, they read names should be globally unique, but
  if they have been stored e.g. in the SRA, the original read names might have
  been lost.
  """

  def __init__(self, iterators):
    # Iterators
    self.iterators = list(iterators)
    # The rest is initialized upon the first call of __next__()
    # Current reads (one per iterator) and last-advanced iterator index.
    self.reads = None
    self.curr = None
    # Current sequence and position
    self.ump = None
    self.seq = None
    self.pos = None

  def __next__(self):
    # I. Initialize iteration if necessary
    if self.reads is None:
      # Draw one element for every iterator
      self.reads = [self.next_or_none(i) for i in self.iterators]
      self.find_leftmost()

    # II. Find the iterator providing the leftmost read.
    # Can use the same iterator as long as it doesn't move to a new position
    # This is the case if it (1) hasn't finished, (2) has neither moved from
    # mapped to unmapped reads, nor to a new chromosome or position.
    r = self.reads[self.curr]
    if r is None or self.ump_seq_pos(r) > (self.ump, self.seq, self.pos):
      self.find_leftmost()
    elif self.ump_seq_pos(r) < (self.ump, self.seq, self.pos):
      raise RuntimeError("input in order of ascending genomic coordinates (%s < %s)" %
                         (self.ump_seq_pos(r), (self.ump, self.seq, self.pos)))

    # III. Return read and advance iterator
    r = self.reads[self.curr]
    self.reads[self.curr] = self.next_or_none(self.iterators[self.curr])
    return r

  @classmethod
  def ump_seq_pos(cls, read):
    return (read.is_unmapped and (not read.is_paired or read.mate_is_unmapped),
            read.reference_name, read.reference_start)

  @classmethod
  def next_or_none(cls, iterator):
    try:
      return iterator.__next__()
    except StopIteration:
      return None

  def find_leftmost(self):
    # Scan current reads and find leftmost
    ump, seq, pos, index = None, None, None, None
    for i, r in enumerate(self.reads):
      if r is not None and (ump is None or self.ump_seq_pos(r) < (ump, seq, pos)):
        ump, seq, pos, index = *self.ump_seq_pos(r), i
    # Update state, raise StopIteration if no more reads
    if index is None:
      self.ump, self.seq, self.pos, self.curr = None, None, None, None
      raise StopIteration()
    else:
      if self.seq is not None and (self.ump, self.seq, self.pos) > (ump, seq, pos):
        raise RuntimeError(
            "input not in order of ascending genomic coordinates")
      self.ump, self.seq, self.pos, self.curr = ump, seq, pos, index


class injectable_iterator:
  """
  Iterator which allows objects to be dynamically injected during the iteration

  At any time, additional objects can be injected with the `inject` member,
  and are then returned by the next call to `__next__()`. After draining
  the list of injected objects, the iteration continues with the next
  object from the input iterator
  """

  def __init__(self, iterator):
    self.iterator = iterator
    self.queue = list()
    
  def __iter__(self):
    return self
    
  def inject(self, object, flag = None):
    self.queue.append((object, flag))
    
  def __next__(self):
    if self.queue:
      return self.queue.pop(0)
    else:
      return (self.iterator.__next__(), None)


class reorder_buffer:
  """
  Stream objects in a particular order even if some objects become available later
  
  Instead of an object, a placeholder can be added to the stream. Further
  objects are then buffered until an actual object is supplied to fill the
  reserved place. At the point, first the fill-in object and then all buffered
  objects are written to the output stream
  """
  
  object_placeholder = namedlist('reorder_buffer_object_placeholder', ['reorder_buffer', 'slot'])
  vacant = object()
  dont_fill = object()

  def __init__(self, writer):
    self.postponed = list()
    self.writer = writer
    
  def write(self, object):
    # If there are postponed objects, enqueue new object at the end
    # Otherwise, output immediately.
    if self.postponed:
      self.postponed.append(object)
    else:
      self.writer.write(object)
  
  def placeholder(self):
    # Create placeholder object and put onto postponed list
    ph = reorder_buffer.object_placeholder(reorder_buffer=self, slot=reorder_buffer.vacant)
    self.postponed.append(ph)
    return ph
  
  def fillin(self, object, placeholder):
    # Validate placeholder
    if placeholder.__class__ != reorder_buffer.object_placeholder:
      raise TypeError('placeholder must be of type reorder_buffer.object_placeholder')
    if id(placeholder.reorder_buffer) != id(self):
      raise ValueError('placeholder belongs to different reorder_buffer')
    if id(placeholder.slot) != id(reorder_buffer.vacant):
      raise ValueError('placeholder was filled in already')
    # Put onto postponed list at correct position
    placeholder.slot = object
    # Write out filled-in slots as far as possible
    self.process_postponed()

  def dontfill(self, placeholder):
    self.fillin(reorder_buffer.dont_fill, placeholder=placeholder)

  def nr_postponed(self):
    return len(self.postponed)

  def process_postponed(self):
    # Write postponed objects up to the first postponed placeholder
    while self.postponed:
      head = self.postponed[0]
      if head.__class__ != reorder_buffer.object_placeholder:
        self.postponed.pop(0)
        self.writer.write(head)
      elif id(head.slot) == id(reorder_buffer.vacant):
        break
      elif id(head.slot) == id(reorder_buffer.dont_fill):
        self.postponed.pop(0)
      else:
        self.postponed.pop(0)
        self.writer.write(head.slot)

  def close(self):
    if self.postponed:
      raise RuntimeError('cannot close reorder_buffer containing placeholders')
    self.writer.close()
