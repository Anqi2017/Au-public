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

  def print_range(self):
    print self.chr+":"+str(self.start)+"-"+str(self.end)

  def overlaps(self,in_genomic_range):
    if self.chr != in_genomic_range.chr:
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
    return False
