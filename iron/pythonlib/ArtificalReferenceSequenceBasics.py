import sys, re, zlib, base64
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
    self.conversion_string = None # Gets set by reading an ARS name or entry, or setting bounds
    self.ars_name = None  # sequence that contains the name and conversion string
    self.name = None # Optionally, A custom name for this ARS (becomes encoded in the ars_name)
    self.bounds = None
    self.sequence = None
  def set_ars_name(self,ars_name):
    self.ars_name = ars_name
    [self.conversion_string, self.name] = decode_ars_name(ars_name)
    self.set_conversion_string(self.conversion_string)
  def set_conversion_string(self,conversion_string):
    self.conversion_string = conversion_string
    self.ars_name = encode_ars_name(conversion_string,self.name)
    self.bounds = []
    for part in conversion_string.split('/'):
      m = re.match('^([^,]+),(\d+)-(\d+)\|([+-])$',part)
      self.bounds.append(Bed(m.group(1),int(m.group(2)),int(m.group(3)),m.group(4)))
  def set_bounds(self,inbounds):
    self.bounds = inbounds
    self.bounds_to_conversion_string() # go ahead and set the name
  #Read in the name and the sequence (with no line breaks or spaces)
    self.ars_name = encode_ars_name(self.conversion_string)
  def read_entry(self,ars_name,sequence):
    [conv_string, entry_name] = decode_ars_name(ars_name)
    self.set_conversion_string(conv_string)
    self.name = entry_name
    self.safe_name = ars_safe_name
    self.sequence = sequence
  def set_name(self, entry_name):
    self.name = entry_name
    if self.conversion_string: # if we have a conversion string already refresh our ars_name
      self.ars_name = encode_ars_name(self.conversion_string,self.name)
  def bounds_to_conversion_string(self):
    self.conversion_string = ''
    for b in self.bounds:
      arr = b.get_bed_array()
      self.conversion_string += arr[0]+','+str(arr[1])+'-'+str(arr[2])+'|'+str(arr[3])+'/'
    self.conversion_string = self.conversion_string[:-1]
    self.ars_name = encode_ars_name(self.conversion_string)
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
    if self.ars_name and self.sequence:
      return '>'+self.ars_name+"\n"+self.sequence+"\n"
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

def encode_ars_name(conversion_string,name=''):
  if not name:
    name = ''
  compressed_string = zlib.compress(conversion_string+"\t"+name,9)
  enc_string = base64.b32encode(compressed_string)
  return 'ARS_'+enc_string.rstrip('=')

def decode_ars_name(safename):
  frag = safename.lstrip('ARS_')
  padding = (8-(len(frag) % 8)) % 8
  c = base64.b32decode(frag+'='*padding)
  [conv,name] = zlib.decompress(c).split("\t")
  if name == '':
    name = None
  return [conv,name]
