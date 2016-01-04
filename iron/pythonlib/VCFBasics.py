import sys, re
from RangeBasics import Bed

class VCF:
  def __init__(self,inline=None):
    self.entry = None
    self.payload = None
    self.range = None
    if inline:
      self.entry = line_to_entry(inline)
      self.range = Bed(self.value('chrom'),self.value('pos')-1,self.value('pos'))
  def set_payload(self,inpay):
    self.payload = [inpay]
  def get_payload(self):
    return self.payload[0]
  def value(self,inkey):
    if inkey not in self.entry:
      sys.stderr.write("ERROR "+inkey+" not in entry\n")
      sys.exit()
    return self.entry[inkey]
  def get_phased_genotype(self):
    if not self.value('remainder'):
      sys.stderr.write("Warning expected something in the remainder\n")
      return False
    m = re.search('([01])\|([01])',self.value('remainder'))
    if not m: return False
    left = int(m.group(1))
    right = int(m.group(2))
    alleles = [self.value('ref'),self.value('alt')]
    return [alleles[left],alleles[right]]
  def is_snp(self):
    if len(self.value('ref'))==1 and len(self.value('alt'))==1:
      return True
    return False

def line_to_entry(line):
  f = line.rstrip().split("\t")
  v = {}
  v['chrom'] = f[0]
  v['pos'] = int(f[1]) # 1-indexed
  v['id'] = f[2]
  v['id_array'] = f[2].split(';')
  v['ref'] = f[3]
  v['alt'] = f[4]
  v['qual'] = float(f[5])
  v['filter'] = f[6]
  v['info'] = f[7]
  v['remainder'] = None
  if len(f) > 8:
    v['remainder'] = "\t".join(f[8:])
  return v

