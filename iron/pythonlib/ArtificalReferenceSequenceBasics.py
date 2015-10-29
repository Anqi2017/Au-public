import sys, re
from SequenceBasics import rc
from RangeBasics import Bed

#This class will support the creation of an artifical reference sequence
#and then the conversion of the coordinates of that sequence

# The ARS name is suffient to handle any desired coordinate conversions
# Getting the actual sequence requires having been provided the sequence
# or the original reference
# Pre: must either set_name or set_bounds or read_entry
class ARS:
  def __init__(self):
    self.name = None
    self.bounds = None
    self.sequence = None
  def set_name(self,ars_name):
    self.name = ars_name
    self.bounds = []
    for part in ars_name.split('/'):
      m = re.match('^([^,]+),(\d+)-(\d+)\|([+-])$',part)
      self.bounds.append(Bed(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4)))
  def set_bounds(self,inbounds):
    self.bounds = inbounds
    self.bounds_to_name() # go ahead and set the name
  #Read in the name and the sequence (with no line breaks or spaces)
  def read_entry(self,ars_name,sequence):
    self.set_name(ars_name)
    self.sequence = sequence
  def bounds_to_name(self):
    self.name = ''
    for b in self.bounds:
      arr = b.get_bed_array()
      self.name += arr[0]+','+str(arr[1])+'-'+str(arr[2])+'|'+str(arr[3])+'/'
    self.name = self.name[:-1]
  # Pre: Requires that bounds already be set
  def set_sequence_from_original_reference_hash(self,ref_hash):
    if self.bounds:
      self.construct_sequences(ref_hash)
    else:
      sys.stderr.write("ERROR: bounds must be set before setting the sequence\n")
      sys.exit()
    return
  def construct_sequences(self,ref_hash):
    self.sequence = ''
    for b in self.bounds:
      arr = b.get_bed_array()
      if re.search(',',arr[0]) or re.search('/',arr[0]) or re.search('\|',arr[0]):
        sys.stderr.write("ERROR: original reference chromosome cannot have a comma, forward slash or vertical bar in its name\n")
        sys.exit()
      if arr[0] not in ref_hash:
        sys.stderr.write("ERROR: sequence "+str(arr[0])+" not found in reference\n")
        sys.exit()
      seq = ''
      if arr[3] == '+':
        seq = ref_hash[arr[0]][arr[1]:arr[2]]
      elif arr[3] == '-':
        seq = rc(ref_hash[arr[0]][arr[1]:arr[2]])
      else:
        sys.stderr.write("ERROR: no direction set in bed\n")
      self.sequence += seq.upper()
    return
  def get_fasta(self):
    if self.name and self.sequence:
      return '>'+self.name+"\n"+self.sequence+"\n"
    else:
      return None

  # Pre:  A 1-indexed ARS_coordiante
  # Post: A 1-indexed genomic coordinate in an array [chr, coord, strand] 
  def convert_ARS_to_genomic_coordinate(self,ARS_coordinate):
    if ARS_coordinate < 1: return None
    covered = 0
    for b in self.bounds:
      blen = b.length()
      if ARS_coordinate <= covered+blen:
        if b.get_direction() == '+':
          num = ARS_coordinate-covered
          return [b.chr,num+b.start-1,b.direction]
        elif b.get_direction() == '-':
          num = ARS_coordinate-covered
          return [b.chr,b.end-(num-1),b.direction]
      covered += blen
    return None

  # Pre:  A 1-indexed genomic coordinate
  # Post: A 1-indexed ARS coordiante [coord,strand]
  def convert_genomic_to_ARS_coordinate(self,genomic_coordinate):
    covered = 0
    for b in self.bounds:
      blen = b.length()
      if genomic_coordinate >= b.start and genomic_coordinate <= b.end:
        if b.direction =='+':
          return [covered+(genomic_coordinate-b.start+1),'+']
        elif b.direction == '-':
          return [covered+(b.end-genomic_coordinate+1),'-']
      covered += blen
    return None

