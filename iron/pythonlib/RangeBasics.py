import sys

class GenomicRangeDictionary:
  def __init__(self):
    self.members = []

  def length(self):
    return len(self.members)

  def get_values(self):
    v = set()
    for m in self.members:
      v.add(m[1])
    return v

  def get_range_list(self):
    v = []
    for m in self.members:
      v.append(m[0])
    return v

  def add(self,genomic_range_key,payload):
    temp = []
    temp.append(genomic_range_key)
    temp.append(payload)
    self.members.append(temp)

  # access based on key
  def get_overlapped(self,genomic_range):
    r = GenomicRangeDictionary()
    for m in self.members:
      if m[0].overlaps(genomic_range):
        r.members.append(m)
    return r

  def remove_overlapped(self,genomic_range):
    r = GenomicRangeDictionary()
    newmembers = []
    for m in self.members:
      if not m[0].overlaps(genomic_range):
        newmembers.append(m)
    self.members = newmembers

class GenomicRange:
  def __init__(self,chr,start,end):
    self.chr = str(chr)
    self.start = int(start)
    self.end = int(end)
   
  def length(self):
    return self.end-self.start+1

  def equals(self,gr):
    if self.chr == gr.chr and self.start == gr.start and self.end == gr.end:
      return True
    return False

  def print_range(self):
    print self.chr+":"+str(self.start)+"-"+str(self.end)

  def overlaps(self,in_genomic_range):
    if self.chr != in_genomic_range.chr:
      return False
    if self.end < in_genomic_range.start:
      return False
    if in_genomic_range.end < self.start:
      return False
    if self.start > in_genomic_range.end:
      return False
    if in_genomic_range.start > self.end:
      return False
    if self.start <= in_genomic_range.start and self.end >= in_genomic_range.start:
      return True
    if self.start <= in_genomic_range.end and self.end >= in_genomic_range.end:
      return True
    if self.start >= in_genomic_range.start and self.end <= in_genomic_range.end:
      return True
    if self.start <= in_genomic_range.start and self.end >= in_genomic_range.end:
      return True
    if in_genomic_range.start <= self.start and in_genomic_range.end >= self.start:
      return True
    if in_genomic_range.start <= self.end and in_genomic_range.end >= self.end:
      return True
    sys.stderr.write("overlaps: unprogrammed error\n")
    return False

  def overlap_size(self,in_genomic_range):
    if self.chr != in_genomic_range.chr:
      return 0
    if self.end < in_genomic_range.start:
      return 0
    if in_genomic_range.end < self.start:
      return 0
    if self.start > in_genomic_range.end:
      return 0
    if self.start >= in_genomic_range.start and self.end <= in_genomic_range.end:
      return self.end-self.start+1
    if self.start <= in_genomic_range.start and self.end >= in_genomic_range.end:
      return in_genomic_range.end-in_genomic_range.start+1
    if self.start <= in_genomic_range.start and self.end >= in_genomic_range.start:
      return self.end-in_genomic_range.start+1
    if self.start <= in_genomic_range.end and self.end >= in_genomic_range.end:
      return in_genomic_range.end-self.start+1
    if in_genomic_range.start <= self.start and in_genomic_range.end >= self.start:
      return in_genomic_range.end-self.start+1
    if in_genomic_range.start <= self.end and in_genomic_range.end >= self.end:
      return self.end-in_genomic_range.start+1
    sys.stderr.write("overlap_size: unprogrammed error\n")
    return 0
