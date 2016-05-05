from Bio.Sequence import rc
from Bio.Simulation.RandomSource import RandomSource

class MakeErrors:
  def __init__(self,rand=None,seed=None):
    if rand:
      self.random = rand
    else:
      self.random = RandomSource()
      if seed: self.random = RandomSource(seed)
    #### context information ####
    self._before_base = None
    self._after_base = None
    #### set the reference base to change for del,mismatch ###
    self._observed_base = None
    #### set waht to change base to for ins or mismatch
    self._modified_base = None

  def set_before_context(self,base):
    self._before_base = base
  def set_after_context(self,base):
    self._after_base = base
  def set_observed_base(self,base):
    self._observed_base = base
  def set_modified_base(self,base):
    self._modified_base = base

  def random_substitution(self,sequence,rate):
    seq = ''
    for i in range(len(sequence)):
      # check context
      prev = None
      if i >= 1: prev = sequence[i-1]
      next = None
      if i < len(sequence)-1: next = sequence[i+1]
      if self._before_base and (not prev or prev != self._before_base): 
        seq+=sequence[i]
        continue
      if self._after_base and (not next or next != self._after_base): 
        seq+=sequence[i]
        continue
      if self._observed_base and (sequence[i] != self._observed_base):
        seq+=sequence[i]
        continue

      rnum = self.random.random()
      if rnum < rate:
        if not self._modified_base:
          seq += self.random.different_random_nt(sequence[i])
        else:
          seq += self._modified_base
      else:
        seq += sequence[i]
    return seq

  def random_deletion(self,sequence,rate):
    seq = ''
    for i in range(len(sequence)):
      # check context
      prev = None
      if i >= 1: prev = sequence[i-1]
      next = None
      if i < len(sequence)-1: next = sequence[i+1]
      if self._before_base and (not prev or prev != self._before_base): 
        seq+=sequence[i]
        continue
      if self._after_base and (not next or next != self._after_base): 
        seq+=sequence[i]
        continue
      if self._observed_base and (sequence[i] != self._observed_base):
        seq+=sequence[i]
        continue

      rnum = self.random.random()
      if rnum >= rate:
        seq += sequence[i]
    return seq

  def random_insertion(self,sequence,rate,max_inserts=1):
    seq = ''
    z = 0
    while self.random.random() < rate and z < max_inserts:
      if self._before_base: break # can't do this one
      if self._after_base:
        if self._after_base != sequence[1]: break
      z += 1
      if self._modified_base:
        seq += self._modified_base
      else:
        seq += self.random.random_nt()
    z = 0
    for i in range(len(sequence)):
      # check context
      prev = sequence[i]
      next = None
      if i < len(sequence)-1: next = sequence[i+1]
      if self._before_base and (not prev or prev != self._before_base): 
        seq+=sequence[i]
        continue
      if self._after_base and (not next or next != self._after_base): 
        seq+=sequence[i]
        continue

      seq += sequence[i]
      while self.random.random() < rate and z < max_inserts:
        z+=1
        if self._modified_base:
          seq += self._modified_base
        else:
          seq += self.random.random_nt()
      z = 0
    return seq

  def random_flip(self,sequence):
    if self.random.random() < 0.5:
      return rc(sequence)
    return sequence

class MakeCuts:
  def __init__(self,rand=None,seed=None):
    if rand:
      self.random = rand
    else:
      self.random = RandomSource()
      if seed: self.random = RandomSource(seed)
    self._gauss_min = None
    self._gauss_mu = None
    self._gauss_sigma = None
    self.set_lr_cuts()

  def get_cut(self,seq):
    l = min(len(seq),max(self._gauss_min,int(self.random.gauss(self._gauss_mu,self._gauss_sigma))))
    leeway = len(seq)-l
    start = self.random.randint(0,leeway)
    return seq[start:start+l]

  def set_lr_cuts(self):
    self._gauss_min = 1000
    self._gauss_mu = 4000
    self._gauss_sigma = 500
  def set_sr_cuts(self):
    self._gauss_min = 150
    self._gauss_mu = 290
    self._gauss_sigma = 290
