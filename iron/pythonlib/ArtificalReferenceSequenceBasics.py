import sys, re, zlib, base64
from SequenceBasics import rc
from RangeBasics import Bed
from PSLBasics import PSL
from MultiplePSLBasics import MultiplePSLAlignments

#This class will support the creation of an artifical reference sequence
#and then the conversion of the coordinates of that sequence

# The ARS name is suffient to handle any desired coordinate conversions
# Getting the actual sequence requires having been provided the sequence
# or the original reference
# Pre: must either set_name or set_bounds or read_entry
class ARS:
  def __init__(self,ref=None,name=None,conversion_string=None):
    self.name = None # Optionally, A custom name for this ARS (becomes encoded in the ars_name)
    self.conversion_string = None # Gets set by reading an ARS name or entry, or setting bounds
    self.ars_name = None  # sequence that contains the name and conversion string
    self.bounds = None
    self.sequence = None
    if name: self.set_name(name)
    if conversion_string: 
      self.set_conversion_string(conversion_string)
    self.ref_hash = ref # hash of reference genome
  def get_ars_name(self):
    return self.ars_name
  def get_conversion_string(self):
    return self.conversion_string
  def set_ref_hash(self,ref):
    self.ref_hash = ref
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
    if self.ref_hash:
      self.set_sequence_from_original_reference_hash()
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
    if self.name: 
      self.ars_name = encode_ars_name(self.conversion_string)
    else:
      self.ars_name = encode_ars_name(self.conversion_string,self.name)
  # Pre: Requires that bounds already be set
  def set_sequence_from_original_reference_hash(self):
    if self.bounds:
      self.construct_sequences(self.ref_hash)
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
    if not self.sequence: self.set_sequence_from_original_reference_hash()
    if self.ars_name and self.sequence:
      return '>'+self.ars_name+"\n"+self.sequence+"\n"
    else:
      return None

  # Pre:  A 1-indexed ARS_coordiante
  # Post: A 1-indexed genomic coordinate in an array [chr, coord, strand] 
  def convert_ARS_to_genomic_coordinate(self,ARS_coordinate,in_strand='+'):
    if ARS_coordinate < 1: return None
    covered = 0
    for b in self.bounds:
      blen = b.length()
      if ARS_coordinate <= covered+blen:
        if b.get_direction() == '+':
          num = ARS_coordinate-covered
          if in_strand == '+':
            return [b.chr,num+b.start-1,b.direction]
          else:
            return [b.chr,num+b.start-1,flipped(b.direction)]
        elif b.get_direction() == '-':
          num = ARS_coordinate-covered
          if in_strand == '+':
            return [b.chr,b.end-(num-1),b.direction]
          else:
            return [b.chr,b.end-(num-1),flipped(b.direction)]
      covered += blen
    return None

  # Pre:  A chromosome
  #       A 1-indexed genomic coordinate
  # Post: A 1-indexed ARS coordiante [coord,strand]
  def convert_genomic_to_ARS_coordinate(self,chrom,genomic_coordinate):
    covered = 0
    for b in self.bounds:
      blen = b.length()
      if b.chr == chrom:
        if genomic_coordinate >= b.start and genomic_coordinate <= b.end:
          if b.direction =='+':
            return [covered+(genomic_coordinate-b.start+1),'+']
          elif b.direction == '-':
            return [covered+(b.end-genomic_coordinate+1),'-']
      covered += blen
    return None
  
  #Convert an ARS alignment psl to a genomic psl
  def convert_ARS_to_genomic_psl(self,psl,maximum_intron=400000):
    newpsl = psl.copy()
    coords = []
    for i in range(0,len(psl.value('blockSizes'))):
      for j in range(0,psl.value('blockSizes')[i]):
        k = psl.value('tStarts')[i]+j
        m = psl.value('qStarts')[i]+j
        v = self.convert_ARS_to_genomic_coordinate(k+1,in_strand=psl.value('strand'))
        if not v: return None
        v.append(m+1)
        v.append(psl.value('strand'))
        coords.append(v)
        #name = self.get_conversion_string() # use this if they forgot to set a name
        #if self.name:  name = self.name
        name = psl.value('qName')
        #print name + "\t" + self.get_conversion_string()+"\t"+str(v)+"\t"+psl.value('strand')
    psl_lines = crush_coords(coords,maximum_intron,name)
    mpsl = MultiplePSLAlignments()
    for psl_line in psl_lines:
      mpsl.add_entry(PSL(psl_line))
    for i in range(0,mpsl.entry_count()): mpsl.entries[i].recalculate_stats()
    return mpsl

def crush_coords(coords,maximum_intron,qname):
  sets = []
  # these all come in ordered by query
  qend = 0
  for coord in coords:
    [chr,target_coord, target_strand, query_coord, query_strand] = coord
    if new_locus(sets,coord,maximum_intron):
      sets.append([])
    sets[-1].append(coord)
    qend = query_coord
  psls = []
  for cset in sets:
    chrom = cset[0][0]
    strand = '+'
    if cset[0][2] != cset[0][4]:  strand = '-'
    if strand == '-': 
      cset.reverse()
      for i in range(0,len(cset)): cset[i][3] = qend-cset[i][3]+1
    blocks = []
    for coord in cset:
      if len(blocks) == 0: blocks.append([])
      if len(blocks[-1])==0: 
        blocks[-1].append(coord)
        continue
      elif blocks[-1][-1][1] != coord[1]-1: # its a new set
        blocks.append([])
      blocks[-1].append(coord)
      continue
    blocksizes =  [x[-1][1]-x[0][1]+1 for x in blocks] 
    targetbounds =  [x[0][1]-1 for x in blocks] 
    querybounds = [x[0][3]-1 for x in blocks] 
    c = ''
    c += str(sum(blocksizes))+"\t"
    c += "0"+"\t"
    c += "0"+"\t"
    c += "0"+"\t"
    c += "0"+"\t"
    c += "0"+"\t"
    c += "0"+"\t"
    c += "0"+"\t"
    c += strand+"\t"
    c += qname+"\t"
    c += str(qend)+"\t"
    c += str(querybounds[0])+"\t"
    c += str(querybounds[-1]+blocksizes[-1])+"\t"
    c += chrom+"\t"
    c += str(targetbounds[-1]+blocksizes[-1])+"\t"
    c += str(targetbounds[0])+"\t"
    c += str(targetbounds[-1]+blocksizes[-1])+"\t"
    c += str(len(blocksizes))+"\t"
    c += ','.join([str(x) for x in blocksizes])+','+"\t"
    c += ','.join([str(x) for x in querybounds])+','+"\t"
    c += ','.join([str(x) for x in targetbounds])+','+"\t"
    psls.append(c)
  return psls

def new_locus(sets,coord,maximum_intron):
  if len(sets) == 0: return True
  if len(sets[-1]) == 0: return False
  if sets[-1][-1][0] != coord[0]: return True
  if sets[-1][-1][2] != coord[2]: return True
  if maximum_intron >= 0:
    if abs(sets[-1][-1][1]-coord[1]) > maximum_intron: return True
  return False

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

def flipped(dir):
  if dir == '+':
    return '-'
  return '+'

class ARS_conversion_string_factory:
  def __init__(self):
    self.bounds = ''
  def add_bounds(self,input_bed,dir=None):
    if input_bed.direction == None and dir == None:
      sys.stderr.write("ERROR direction must be set for range\n")
      sys.exit()
    if input_bed.direction and dir:
      if input_bed.direction != dir: 
        sys.stderr.write("ERROR,direction mismatch?\n")
    if input_bed.direction: dir = input_bed.direction
    if len(self.bounds) > 0: 
      self.bounds += '/'
    self.bounds += input_bed.chr+','+str(input_bed.start)+'-'+str(input_bed.end)+'|'+dir
    return
  def set_conversion_string(self,inbounds): self.bounds = inbounds
  def get_conversion_string(self): return self.bounds
