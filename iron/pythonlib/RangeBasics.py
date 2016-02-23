import sys, re
# These classes are to help deal with genomic coordinates and 
# this associated with those coordinates.

class GenomicRangeDictionary:
  def __init__(self):
    self.members = []

  def length(self):
    return len(self.members)

  def get_range_list(self):
    v = []
    for m in self.members:
      v.append(m)
    return v

  def add(self,genomic_range):
    self.members.append(genomic_range)

  # access based on key
  def get_overlapped(self,genomic_range):
    r = GenomicRangeDictionary()
    for m in self.members:
      if m.overlaps(genomic_range):
        r.add(m)
    return r

  def remove_overlapped(self,genomic_range):
    r = GenomicRangeDictionary()
    newmembers = []
    for m in self.members:
      if not m.overlaps(genomic_range):
        newmembers.append(m)
    self.members = newmembers

#These are 1-index for both start and end
class GenomicRange:
  def __init__(self,chr,start,end,dir=None):
    self.chr = str(chr)
    self.start = int(start)
    self.end = int(end)
    self.direction = dir
    self.payload = [] # should be a reference since its an array

  # Copy with the exception of payload.  Thats still a link
  def copy(self):
    n = GenomicRange(self.chr,self.start,self.end,self.direction)
    n.payload = []
    for p in self.payload:
      n.payload.append(p)
    return n
  def get_payload(self):
    return self.payload[0]
  def set_payload(self,inpay):
    if len(self.payload)==0:
      self.payload = [inpay]
      return
    self.payload[0] = inpay
  def get_direction(self):
    return self.direction
  def set_direction(self,dir):
    self.direction = dir

  def length(self):
    return self.end-self.start+1

  def equals(self,gr):
    if self.chr == gr.chr and self.start == gr.start and self.end == gr.end:
      return True
    return False

  def get_range_string(self):
    return self.chr+":"+str(self.start)+"-"+str(self.end)

  def print_range(self):
    print self.get_range_string()
  # These are the 0-indexed start, 1-indexted stop coordinates
  def get_bed_coordinates(self):
    return [self.chr,self.start-1,self.end]

  # These are the 1-indexed coordiantes
  def get_genomic_coordinates(self):
    return [self.chr,self.start,self.end]

  def overlaps(self,in_genomic_range,use_direction=False):
    if self.chr != in_genomic_range.chr:
      return False
    if self.direction != in_genomic_range.direction and use_direction :
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

  def overlaps_with_padding(self,in_genomic_range,padding):
    in_range_padded = GenomicRange(in_genomic_range.chr,max([1,in_genomic_range.start-padding]),in_genomic_range.end+padding)
    return self.overlaps(in_range_padded)

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

  def merge(self,range2,use_direction=False): #merge this bed with another bed
    if self.chr != range2.chr:
      return None
    if self.direction != range2.direction and use_direction:
      return None
    o = GenomicRange(self.chr,min(self.start,range2.start),max(self.end,range2.end),self.direction)
    if use_direction==False: o.direction = None
    return o

  # return 1 if greater than range2
  # return -1 if less than range2
  # return 0 if overlapped
  def cmp(self,range2,overlap_size=0):
    if self.overlaps_with_padding(range2,overlap_size): return 0
    if self.chr < range2.chr: return -1
    elif self.chr > range2.chr: return 1
    if self.end < range2.start: return -1
    elif self.start > range2.end: return 1
    sys.stderr.write("ERROR: cmp function unexpcted state\n")
    sys.exit()
    return 0
  
  #pre another rnage
  #post a list of ranges after removing range2
  #no garuntees on payload
  def subtract(self,range2,use_direction=False):
    outranges = []
    if self.chr != range2.chr:
      outranges.append(self.copy())
      return outranges
    if self.direction != range2.direction and use_direction:
      outranges.append(self.copy())
      return outranges
    if not self.overlaps(range2,use_direction=use_direction):
      outranges.append(self.copy())
      return outranges
    #print self.get_range_string()
    #print range2.get_range_string()
    #print '---'
    if range2.start <= self.start and range2.end >= self.end:
      return outranges #delete all
    if range2.start > self.start: #left side
      nrng = Bed(self.chr,self.start-1,range2.start-1,self.direction)
      outranges.append(nrng)
    if range2.end < self.end: #right side
      nrng = Bed(self.chr,range2.end,self.end,self.direction)
      outranges.append(nrng)
    return outranges

# Pre: Inherits all methods of GenomicRange but modifies the class to use the 0-based start 1-based end style of a bed file
# Essentially, a Bed is just another way of defining a GenomicRange.
class Bed(GenomicRange):
  # Takes as an input the chromosome, 0-based start, 1-based end of a range
  def __init__(self,chrom,start,finish,dir=None):
    self.start = int(start)+1
    self.end = int(finish)
    self.chr = str(chrom)
    self.payload = []
    self.direction = dir
  def copy(self):
    n = Bed(self.chr,self.start-1,self.end,self.direction)
    n.payload = []
    for p in self.payload:
      n.payload.append(p)
    return n
  def get_bed_array(self):
    arr = [self.chr,self.start-1,self.end]
    if self.direction:
      arr.append(self.direction)
    return arr 

class Locus:
  # A Locus is a colloction of GenomicRanges that fall within some distance of one another
  def __init__(self):
    self.range = None
    self.members = []
    self.use_direction = False #If we want to use direction, input ranges must have the direction set
  # Set to true if you want all locus members to share the same direction
  def set_use_direction(self,inbool):
    self.use_direction = inbool
  def add_member(self,grange):
    if self.use_direction and not grange.direction:
      sys.stderr.write("ERROR if using direction then direction of input members must be set\n")
      sys.exit()
    # Get range set properly
    if not self.range:
      self.range = GenomicRange(grange.chr,grange.start,grange.end)
      if self.use_direction:
        self.range.set_direction(grange.get_direction())
    elif self.range.chr != grange.chr:
      sys.stderr.write("WARNING cannot add member with chromosomes are not equal\n")
      return False
    elif self.use_direction and self.range.direction != grange.direction:
      sys.stderr.write("WARNING cannot add member with different directions\n")
      return False
    else:
      if grange.start < self.range.start:  self.range.start = grange.start
      if grange.end > self.range.end: self.range.end = grange.end
    self.members.append(grange)

class Loci:
  # multiple locus combined together when new members are added
  # based on parameters
  def __init__(self):
    self.loci = []
    self.overhang = 0
    self.use_direction = False
    self.verbose = False
    return
  def set_minimum_distance(self,over):
    self.overhang = over
  # Do we want to only combine loci when they have the same direction, if so, set to True
  def set_use_direction(self,inbool):
    self.use_direction = inbool
  # Adds a locus to our loci, but does not go through an update our locus sets yet
  def add_locus(self,inlocus):
    if self.use_direction == True and inlocus.use_direction == False:
      sys.stderr.write("ERROR if using the direction in Loci, then every locus added needs use_direction to be True\n")
      sys.exit()
    self.loci.append(inlocus)
    return
  # Goes through and combines loci until we have one set meeting our overlap definition
  def update_loci(self):
    # Create sub-loci for each chromosome
    lbc = {}
    chroms = sorted([x.range.chr for x in self.loci])
    for chrom in chroms: lbc[chrom] = Loci()
    for x in self.loci: lbc[x.range.chr].add_locus(x)
    for chrom in sorted(lbc.keys()):
      if self.verbose: 
        lbc[chrom].verbose = True
        sys.stderr.write(chrom+"\n")
      lbc[chrom].overhang = self.overhang
      lbc[chrom].use_direction = self.use_direction
      lbc[chrom].merge_down_loci()
    self.loci = []
    for chrom in sorted(lbc.keys()):
      for locus in lbc[chrom].loci:  self.loci.append(locus)
  def merge_down_loci(self):
    old_locus_size = -1
    z = 0
    while len(self.loci) != old_locus_size:
      z+=1
      old_locus_size = len(self.loci)
      locus_size = len(self.loci)
      if self.verbose:
        sys.stderr.write(str(locus_size)+" Combining down loci step "+str(z)+"       \r")
      combined = set()
      for i in range(0,locus_size):
        if i in combined: continue
        for j in range(i+1,locus_size):
          if self.loci[i].range.overlaps_with_padding(self.loci[j].range,self.overhang):
            if self.use_direction and self.loci[i].range.direction != self.loci[j].range.direction:  continue
            for obj in self.loci[j].members:
              self.loci[i].add_member(obj)
            combined.add(j)
            break
      newloci = []
      for i in range(0,locus_size):
        if i not in combined:
          newloci.append(self.loci[i])
      self.loci = newloci
    if self.verbose:
      sys.stderr.write("Finished combining down "+str(len(self.loci))+" loci in "+str(z)+" steps   \n")
    return

#pre an array of ranges
#post a sorted array of ranges
def sort_ranges(inranges):
  if not inranges: return
  outranges = []
  if len(inranges)==0: return outranges
  v = {}
  for rng in inranges:
    if rng.chr not in v: v[rng.chr] = {}
    if rng.start not in v[rng.chr]: v[rng.chr][rng.start] = {}
    if rng.end not in v[rng.chr][rng.start]: v[rng.chr][rng.start][rng.end] = []
    v[rng.chr][rng.start][rng.end].append(rng)
  for chr in sorted(v.keys()):
    for start in sorted(v[chr].keys()):
      for end in sorted(v[chr][start].keys()):
        for e in v[chr][start][end]:
          outranges.append(e)
  return outranges
  
#Pre: list of bed tools, whether or not they are already sorted
#Post: flattend range list of ranges where if they overlapped, they are now joined
#      (not yet) The new range payloads will be the previous ranges
def merge_ranges(inranges,already_sorted=False):
  if not already_sorted: inranges = sort_ranges(inranges)
  prev = None
  outputs = []
  merged = False
  for rng in inranges:
    nrng = rng.copy()
    nrng.set_payload([])
    nrng.get_payload().append(rng)
    merged = False
    if prev:
      if rng.overlaps(prev):
        nrng = nrng.merge(prev,use_direction=False)
        nrng.set_payload(prev.get_payload())
        nrng.get_payload().append(rng)
        merged = True
      else:
        outputs.append(prev)
    prev = nrng
  if not merged: outputs.append(prev)
  return sort_ranges(outputs)

def pad_ranges(inranges,padding,chr_ranges=None):
  if not inranges: return
  outranges = []
  if len(inranges) == 0: return outranges
  chr = {}
  if chr_ranges:
    for b in chr_ranges:
      chr[b.chr] = b
  for rng in inranges:
    newstart = rng.start - padding
    newend = rng.end + padding
    if rng.chr in chr:
      if newstart < chr[rng.chr].start: newstart = chr[rng.chr].start
      if newend > chr[rng.chr].end: endstart = chr[rng.chr].end
    nrng = rng.copy()
    nrng.start = newstart
    nrng.end = newend
    outranges.append(nrng)
  return sort_ranges(outranges)

def subtract_ranges(r1s,r2s,already_sorted=False):
  if not already_sorted:
    r1s = sort_beds(r1s)
    r2s = sort_beds(r2s)
  left = r1s[:]
  right = r2s[:]
  curleft = None
  curright = None
  outputs = []
  while len(left) > 0 and len(right) > 0:
    if len(right) == 0:  outputs.append(left.pop(0))
    if len(left) == 0: right.pop(0)
    c = left[0].cmp(right[0])
    if c == 0:
      s = left[0].subtract(right[0])
      left.pop(0)
      right.pop(0)
      #print '---'
      for i in reversed(s):
        #print 'insert'  
        left.insert(0,i)
    elif c < 0:
      #left is altogether smaller
      outputs.append(left.pop(0))
    elif c > 0:
      #right is altogether smalelr
      right.pop(0)
  return sort_ranges(outputs)

def string_to_genomic_range(rstring):
  m = re.match('([^:]+):(\d+)-(\d+)',rstring)
  if not m: 
    sys.stderr.write("ERROR: problem with range string "+rstring+"\n")
  return GenomicRange(m.group(1),int(m.group(2)),int(m.group(3)))
