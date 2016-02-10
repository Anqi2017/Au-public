from SequenceBasics import rc
# Class for describing the Kmer properties of sequences

class Kmer:
  def __init__(self,ksize=10):
    #self.sequence = None
    self.payload = None # Allow some information in list form associated with the Kmer
    self.ksize = ksize
    self.mers = {}
    return
  def add_sequence(self,seq):
    seq = seq.upper()
    for i in range(len(seq)+1-self.ksize):
      mer = seq[i:i+self.ksize]
      if mer not in self.mers:
        self.mers[mer] = 0
      self.mers[mer] += 1
      if rc(mer) not in self.mers:
        self.mers[rc(mer)] = 0
      self.mers[rc(mer)] += 1
    return
  def load_kmer(self,inkmers):
    self.ksize = inkmers.ksize
    for i in inkmers.mers:
      self.mers[i] = inkmers.mers[i]
    return
  # Take another Kmer object and return the kmer sequences that are in both
  def overlaps(self,inkmer):
    ov = set()
    for i in inkmer.mers:
      if i in self.mers:
        ov.add(i)
    #if len(ov) == 0: return None
    return ov
  def add_kmer(self,inkmer):
    for i in inkmer.mers:
      if i not in self.mers:
        self.mers[i] = inkmer.mers[i]
      else: 
        self.mers[i]+= inkmer.mers[i]
  def set_payload(self,inpay):
    self.payload = [inpay]
  def get_payload(self):
    return self.payload[0]

def zero_count_Kmer(size):
  k = Kmer(size)
  cs = ['A','C','T','G']
  seqs = cs[:]
  for i in range(1,size):
    newseqs = []
    for j in range(0,len(seqs)):
      cur = seqs[j]
      for c in cs:
        newseqs.append(seqs[j]+c)
    seqs = newseqs
  for seq in seqs:
    k.mers[seq] = 0
  return k

class MemoryKmer():
  def __init__(self,size):  
    self.mers = {}
    self.values = []
    self.ksize = size
    cs = ['A','C','T','G']
    seqs = cs[:]
    for i in range(1,size):
      newseqs = []
      for j in range(0,len(seqs)):
        cur = seqs[j]
        for c in cs:
          newseqs.append(seqs[j]+c)
      seqs = newseqs
    z = 0
    for seq in seqs:
      self.mers[seq] = z
      z += 1
      self.values.append(0)
    return
