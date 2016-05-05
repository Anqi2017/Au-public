import sys, random, string
from Bio.Range import GenomicRange
from Bio.Sequence import rc

class Transcript:
  def __init__(self):
    self.exons = []
    self.junctions = []
    self._direction = None
    self._transcript_name = None
    self._gene_name = None

  def get_length(self):
    return sum([x.get_length() for x in self.exons])

  def set_strand(self,dir):
    self._direction = dir
  def get_strand(self):
    return self._direction

  #greedy return the first chromosome in exon array
  def get_chrom(self):
    if len(self.exons)==0: 
      sys.stderr.write("WARNING can't return chromsome with nothing here\n")
      return None
    return self.exons[0].get_range().chr

  # Pre: A strcutre is defined
  #      The Sequence from the reference
  def get_sequence(self,ref_dict):
    strand = '+'
    if not self._direction:
      sys.stderr.write("WARNING: no strand information for the transcript\n")
    if self._direction: strand = self._direction
    chr = self.get_chrom()
    seq = ''
    for e in [x.get_range() for x in self.exons]:
      seq += ref_dict[chr][e.start-1:e.end]
    if strand == '-':  seq = rc(seq)
    return seq.upper()

  def get_gpd_line(self,transcript_name=None,gene_name=None,strand=None):
    tname = self._transcript_name
    gname = self._gene_name
    dir = self._direction
    if not tname: tname = transcript_name
    if not gname: gname = gene_name
    if not dir: dir = strand
    if not tname or not gname or strand:
      sys.stderr.write("ERROR:  transcript name and gene name and direction must be set to output a gpd line or use get_fake_gpd_line()\n")
    out = ''
    out += tname + "\t"
    out += gname + "\t"
    out += self.exons[0].rng.chr + "\t"
    out += dir + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(len(self.exons)) + "\t"
    out += str(','.join([str(x.rng.start-1) for x in self.exons]))+','+"\t"
    out += str(','.join([str(x.rng.end) for x in self.exons]))+','
    return out

  def set_gene_name(self,name):
    self._gene_name = name
  def get_gene_name(self):
    return self._gene_name
  def set_transcript_name(self,name):
    self._transcript_name = name
  def get_transcript_name(self):
    return self._transcript_name

  def get_fake_gpd_line(self):
    rlen = 8
    name = ''.join(random.choice(string.letters+string.digits) for i in range(0,rlen))
    out = ''
    out += name + "\t"
    out += name + "\t"
    out += self.exons[0].rng.chr + "\t"
    out += '+' + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(self.exons[0].rng.start-1) + "\t"
    out += str(self.exons[-1].rng.end) + "\t"
    out += str(len(self.exons)) + "\t"
    out += str(','.join([str(x.rng.start-1) for x in self.exons]))+','+"\t"
    out += str(','.join([str(x.rng.end) for x in self.exons]))+','
    return out

  def get_junctions_string(self):
    return ';'.join([x.get_range_string() for x in self.junctions])

  def junction_overlap(self,tx,tolerance=0):
    return Transcript.JunctionOverlap(self,tx,tolerance)

  class JunctionOverlap:
    def __init__(self1,tx_obj1,tx_obj2,tolerance=0):
      self1.tx_obj1 = tx_obj1
      self1.tx_obj2 = tx_obj2
      self1.tolerance = tolerance
      self1.overs = [] # gets set by calculate_overlap()
      self1.calculate_overlap()
      if len(self1.overs) == 0: return # nothing to analyze
      self1.analyze_overs()

    def __nonzero__(self1):
      if len(self1.overs) > 0: return True
      return False

    def match_junction_count(self1):
      return len(self1.overs)

    # Return value if tx_obj2 is a complete subset of tx_obj1 or tx_obj1 is a complete subset of tx_obj2
    # Return 1: Full overlap (mutual subests) 
    # Return 2: two is a subset of one
    # Return 3: one is a subset of two
    # Return False if neither is a subset of the other
    def is_subset(self1):
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0: # make sure they are consecutive if more than one
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      onecov = self1.start1 and self1.end1
      twocov = self1.start2 and self1.end2
      if onecov and twocov:
        return 1
      elif twocov: return 2
      elif onecov: return 3
      return False

    def is_full_overlap(self1):
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      if self1.start1 and self1.end1 and self1.start2 and self1.end2:
        return True
      return False

    # Return True if the transcripts can be combined together
    def is_compatible(self1):
      if len(self1.overs) == 0: return False
      if len(self1.dif1) > 0:
        if max(self1.dif1) != 1 or max(self1.dif2) != 1: return False
      # If we are still here it is a single run
      if (self1.start1 or self1.start2) and (self1.end1 or self1.end2):
        return True
      return False

    def analyze_overs(self1):
      #check for full overlap first
      self1.dif1 = [self1.overs[i][0]-self1.overs[i-1][0] for i in range(1,len(self1.overs))]
      self1.dif2 = [self1.overs[i][1]-self1.overs[i-1][1] for i in range(1,len(self1.overs))]
      #see if it starts and ends on first or last junction
      self1.start1 = self1.overs[0][0] == 0
      self1.start2 = self1.overs[0][1] == 0
      self1.end1 = self1.overs[-1][0] == len(self1.tx_obj1.junctions)-1
      self1.end2 = self1.overs[-1][1] == len(self1.tx_obj2.junctions)-1
      return

    #Create the array that describes how junctions overlap
    def calculate_overlap(self1):
      overs = []
      for i in range(0,len(self1.tx_obj1.junctions)):
        for j in range(0,len(self1.tx_obj2.junctions)):
          if self1.tx_obj1.junctions[i].overlaps(self1.tx_obj2.junctions[j],self1.tolerance):
            overs.append([i,j])
      self1.overs = overs

class Junction:
  def __init__(self,rng_left=None,rng_right=None):
    self.left = rng_left
    self.right = rng_right
    self.left_exon = None
    self.right_exon = None
  def get_left_exon(self):
    return self.left_exon
  def get_right_exon(self):
    return self.right_exon
  def get_range_string(self):
    return self.left.chr+":"+str(self.left.end)+'/'+self.right.chr+":"+str(self.right.start)
  def set_left(self,rng):
    self.left = rng
  def set_right(self,rng):
    self.right = rng
  #test equality with another junction
  def equals(self,junc):
    if self.left.equals(junc.left): return False
    if self.right.equals(junc.right): return False
    return True
  # see if junction overlaps with tolerance    
  def overlaps(self,junc,tolerance=0):
    if not self.left.overlaps_with_padding(junc.left,tolerance): return False
    if not self.right.overlaps_with_padding(junc.right,tolerance): return False
    return True
  #output 
  # -1 if junc comes before self
  # 1 if junc comes after self
  # 0 if overlaps
  # 2 if else
  def cmp(self,junc,tolerance=0):
    if self.overlaps(junc,tolerance):
      return 0 #equal
    if self.left.chr == junc.right.chr:
      if self.left.start > junc.right.start:
        return -1 #comes before
    if self.right.chr == junc.left.chr:
      if self.right.start < junc.right.start:
        return 1 #comes after
    return 2

  def set_exon_left(self,ex):
    self.left_exon = ex
    ex.right_junc = self
  def set_exon_right(self,ex):
    self.right_exon = ex
    ex.left_junc = self

class Exon:
  def __init__(self,rng=None):
    self.rng = rng
    self.left_junc = None
    self.right_junc = None
    self.start = False
    self.end = False
  def get_range(self):
    return self.rng
  def get_length(self):
    return self.rng.length()
  def set_left_junc(self,junc):
    self.left_junc = junc
    junc.set_right_exon = self
  def set_right_junc(self,junc):
    self.right_junc = junc
    junc.set_left_exon = self
  def set_start(self,boo=True): self.start = boo
  def set_end(self,boo=True): self.end = boo 

# A transcript group is like the fuzzy gpd class we had before
class TranscriptGroup:
  def __init__(self):
    self.junction_groups = [] # These will be more fuzzy defitions
    #self.exons = [] # These will be based on the junctions and individual starts
    self.transcripts = [] # These are the individual transcripts that make up this group

  # Return a representative transcript object
  def get_transcript(self,exon_bounds='max'):
    out = Transcript()
    out.junctions = [x.get_junction() for x in self.junction_groups]
    # get internal exons
    self.exons = []
    for i in range(0,len(self.junction_groups)-1):
      j1 = self.junction_groups[i].get_junction()
      j2 = self.junction_groups[i+1].get_junction()
      e = Exon(GenomicRange(j1.right.chr,j1.right.end,j2.left.start))
      e.set_left_junc(j1)
      e.set_right_junc(j2)
      #print str(i)+" to "+str(i+1)
      out.exons.append(e)
    # get left exon
    left_exons = [y for y in [self.transcripts[e[0]].junctions[e[1]].get_left_exon() for e in self.junction_groups[0].evidence] if y]
    if len(left_exons) == 0:
      sys.stderr.write("ERROR no left exon\n")
      sys.exit()
    e_left = Exon(GenomicRange(out.junctions[0].left.chr,\
                               min([x.get_range().start for x in left_exons]),
                               out.junctions[0].left.start))
    e_left.set_right_junc(out.junctions[0])
    out.exons.insert(0,e_left)
    # get right exon
    right_exons = [y for y in [self.transcripts[e[0]].junctions[e[1]].get_right_exon() for e in self.junction_groups[-1].evidence] if y]
    if len(right_exons) == 0:
      sys.stderr.write("ERROR no right exon\n")
      sys.exit()
    e_right = Exon(GenomicRange(out.junctions[-1].right.chr,\
                               out.junctions[-1].right.end,\
                               max([x.get_range().end for x in right_exons])))
    e_right.set_left_junc(out.junctions[-1])
    out.exons.append(e_right)
    return out

  def add_transcript(self,tx,juntol=0):
    # check existing transcripts for compatability
    for t in self.transcripts:
      ov = t.junction_overlap(tx)
      if ov:
        if not ov.is_compatible: return False
      else: return False # if its not overlapped we also can't add
    self.transcripts.append(tx)
    curr_tx = len(self.transcripts)-1
    #print curr_tx
    # see if there is no junctions yet
    if len(self.junction_groups) == 0:
      for i in range(0,len(tx.junctions)):
        jg = TranscriptGroup.JunctionGroup(self)
        jg.add_junction(curr_tx,i,tolerance=juntol)
        self.junction_groups.append(jg)
    else: # there is already a transcript(s) here to work around
      before = []
      middle = []
      after = []
      for j in range(0,len(tx.junctions)):
        # see if its before the existing set
        cmp = self.junction_groups[0].get_junction().cmp(tx.junctions[j])
        if cmp == -1: before.append(j)
        # see if it goes in the existing set
        for k in range(0,len(self.junction_groups)):
          ov = self.junction_groups[k].get_junction().overlaps(tx.junctions[j],tolerance=juntol) #may need to add a tolerance
          if ov: middle.append([j,k])
        # see if it goes after this set
        cmp = self.junction_groups[-1].get_junction().cmp(tx.junctions[j])
        if cmp == 1: after.append(j)
      # add to the middle values before we disrupt indexing
      #print '---'
      #print len(before)
      #print len(middle)
      #print len(after)
      #print '---'
      for v in middle:
        self.junction_groups[v[1]].add_junction(curr_tx,v[0],tolerance=juntol) 
      #add to the beginning and then the end
      for i in reversed(before):
         jg = TranscriptGroup.JunctionGroup(self)
         jg.add_junction(curr_tx,i,tolerance=juntol)
         self.junction_groups.insert(0,jg)        
      for i in after:
         jg = TranscriptGroup.JunctionGroup(self)
         jg.add_junction(curr_tx,i,tolerance=juntol)
         self.junction_groups.append(jg)        
      #if len(tx.junctions)==0:
      #  jg = TranscriptGroup.JunctionGroup(self)
      #  jg.add_junction(curr_tx,i)
      #  self.junctions.append(jg)

  class JunctionGroup:
    def __init__(self1,outer):
      self1.outer = outer
      self1.evidence = [] # array of evidence that is the 
                         # outer.transcript index
                         # outer.trascript.junction index
      self1.representative_junction = None #calculated as needed
    def get_junction(self1): # return the consensus junction
      if self1.representative_junction:
        return self1.representative_junction
      left_rngs = []
      right_rngs = []
      for j in [self1.outer.transcripts[x[0]].junctions[x[1]] for x in self1.evidence]:
        left_rngs.append(j.left)
        right_rngs.append(j.right)
      left = _mode([x.end for x in left_rngs])
      right = _mode([x.start for x in right_rngs])
      outj = Junction(GenomicRange(left_rngs[0].chr,left,left),GenomicRange(right_rngs[0].chr,right,right))
      self1.representative_junction = outj
      return outj
    def add_junction(self1,tx_index,junc_index,tolerance=0):
      self1.representative_junction = None
      if len(self1.evidence)==0:
        # go ahead and add it
        #j = self1.outer.transcripts[tx_index].junctions[junc_index]
        self1.evidence.append([tx_index,junc_index])
      else:
        # check it and add it
        if not self1.get_junction().overlaps(self1.outer.transcripts[tx_index].junctions[junc_index],tolerance=tolerance):
          sys.stderr.write("WARNING Unable to add junction JunctionGroup\n"+self1.get_junction().get_range_string()+"\n"+self1.outer.transcripts[tx_index].junctions[junc_index].get_range_string()+"\n")
          return False
        self1.evidence.append([tx_index,junc_index])

def _mode(mylist):
  counts = [mylist.count(x) for x in mylist]
  maxcount = max(counts)
  avg = sum([float(x) for x in mylist])/len(mylist)
  #print counts 
  dist = [abs(float(x)-avg) for x in mylist]
  best_list = []
  best_dist = []
  for i in range(0,len(mylist)):
    counts[i] == maxcount
    best_list.append(mylist[i])
    best_dist.append(dist[i])
  abs_best_dist = min(best_dist)
  for i in range(0,len(best_dist)):
    if best_dist[i] == abs_best_dist: 
      return best_list[i]
  sys.stderr.write("Warning: trouble finding best\n")
  return best_list[0]

class Transcriptome:
  def __init__(self,gpd_file=None):
    self.transcripts = []
    if gpd_file:
      from Bio.Format.GPD import GPD
      with open(gpd_file) as inf:
        for line in inf:
          self.transcripts.append(GPD(line))
  def get_transcripts(self):
    return self.transcripts
      
  def add_transcript(self,transcript):
    self.transcripts.append(transcript)

  def __str__(self):
    ostr = ''
    ostr += "Transcriptome containing "+str(len(self.transcripts))+" transcripts "
    ostr += "covering "+str(sum([x.get_length() for x in self.transcripts]))+" bases"
    return ostr
