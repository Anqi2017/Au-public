import random

# emits a random nucleotide sequence
# user can set parameters one at a time as desired
#
# - set parameters
# set_gc_content(float default 0.5)
# set_sequence_length_range(min_len,max_len)
# set_sequence_length(length)
#
# emit()  - outputs the sequence
#
class RandomSequenceEmitter:
  def __init__(self):
    # some random settings
    self.sequence_length_min = 200
    self.sequence_length_max = 200
    self.gc_content = 0.5

  def emit(self):
    thislen = random.randint(self.sequence_length_min,self.sequence_length_max)
    o = ''
    for i in range(0,thislen): o += random_nucleotide_gc(self.gc_content)
    return o

  def set_length_range(self,minl,maxl):
    self.sequence_length_min = minl
    self.sequence_length_max = maxl

  def set_length(self,onelen):
    self.sequence_length_min = onelen
    self.sequence_length_max = onelen

  def set_gc_content(self,gc_content):
    self.gc_content = gc_content

def random_nucleotide_gc(gc_content):
  if random.random() < float(gc_content):
    if random.random() < 0.5: return 'G'
    return 'C'
  else:
    if random.random() < 0.5: return 'A'
    return 'T'
