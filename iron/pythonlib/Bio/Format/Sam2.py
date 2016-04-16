import struct, zlib, sys, re
from cStringIO import StringIO
from string import maketrans
_bam_ops = maketrans('012345678','MIDNSHP=X')
_bam_char = maketrans('abcdefghijklmnop','=ACMGRSVTWYHKDBN')
_bam_value_type = {'c':[1,'<b'],'C':[1,'<B'],'s':[2,'<h'],'S':[2,'<H'],'i':[4,'<i'],'I':[4,'<I']}

# A sam entry
class SAM:
  def __init__(self,line=None,dict=None):
    self._entries = {}
    self._line = None
    return
  #Bam files need a specific override to get_tags and get_cigar that would break other parts of the class if we 
  # access the variables other ways
  #Force tags and cigars to be hidden so we don't accidently change them.
  class PrivateValues:
    def __init__(self):
      self.__tags = None
      self.__cigar = None
    def set_tags(self,tags): self.__tags=tags
    def get_tags(self): return self.__tags
    def set_cigar(self,cigar): self.__cigar=cigar
    def get_cigar(self): return self.__cigar
  def get_tags(self): 
    return self._private_values.get_tags()
  def get_cigar(self): 
    return self._private_values.get_tags()
  def __str__(self):
    return self._line
  def value(self,key):
    return self._entries[key]

# Very much like a sam entry but optimized for access from a bam
# Slows down for accessing things that need more decoding like
# sequence, quality, cigar string, and tags
class BAM(SAM):
  def __init__(self,bin_data,ref_names,fileName=None,blockStart=None,innerStart=None):
    part_dict = _parse_bam_data_block(bin_data,ref_names)
    self._entries = part_dict
    self._line = None
    self._file_position = {'fileName':fileName,'blockStart':blockStart,'innerStart':innerStart} # The most special information about the bam
    self._private_values = BAM.PrivateValues() # keep from accidently accessing some variables other than by methods
    return

  def get_filename(self):
    return self._file_position['fileName']
  def get_block_start(self):
    return self._file_position['blockStart']
  def get_inner_start(self):
    return self._file_position['innerStart']
  def get_file_position_string(self):
    return 'fileName: '+self._file_position['fileName']+" "\
           'blockStart: '+str(self._file_position['blockStart'])+" "\
           'innerStart: '+str(self._file_position['innerStart'])
  def get_tags(self): 
    cur = self._private_values.get_tags()
    if not cur:
      v1,v2 = _bin_to_extra(self._entries['extra_bytes'])
      self._private_values.set_tags(v1) #keep the cigar array in a special palce
      self._entries['remainder'] = v2
    return self._private_values.get_tags()
  def get_cigar(self): 
    cur = self._private_values.get_cigar()
    if not cur:
      v1,v2 = _bin_to_cigar(self._entries['cigar_bytes'])
      self._private_values.set_cigar(v1) #keep the cigar array in a special palce
      self._entries['cigar'] = v2
    return self._private_values.get_cigar()

  def __str__(self):
    if not self._line:
      if self.value('rnext') == self.value('rname'): self._entries['rnext'] = '='
      self._line = self.value('qname')+"\t"+str(self.value('flag'))+"\t"+self.value('rname')+"\t"+str(self.value('pos'))+"\t"+str(self.value('mapq'))+"\t"+self.value('cigar')+"\t"+self.value('rnext')+"\t"+str(self.value('pnext'))+"\t"+str(self.value('tlen'))+"\t"+self.value('seq')+"\t"+self.value('qual')
      if self.value('remainder'):
        self._line += "\t"+self.value('remainder')
    return self._line


  def value(self,key):
    if key not in self._entries:
      if key == 'seq':
        v = _bin_to_seq(self._entries['seq_bytes'])
        self._entries['seq'] = v
        return v
      elif key == 'qual':
        v = _bin_to_qual(self._entries['qual_bytes'])
        self._entries['qual'] = v
        return v
      elif key == 'cigar':
        v1,v2 = _bin_to_cigar(self._entries['cigar_bytes'])
        self._private_values.set_cigar(v1) #keep the cigar array in a special palce
        self._entries['cigar'] = v2
        return self._entries[key]
      elif key == 'remainder':
        v1,v2 = _bin_to_extra(self._entries['extra_bytes'])
        self._private_values.set_tags(v1) #keep the cigar array in a special palce
        self._entries['remainder'] = v2
        return self._entries[key]
    return self._entries[key]

class BAMFile:
  def __init__(self,filename,blockStart=None,innerStart=None,cnt=None):
    self.path = filename
    self.fh = BGZF(filename)
    # start reading the bam file
    self.header_text = None
    self.n_ref = None
    self._read_top_header()
    self.ref_names = []
    self.ref_lengths = {}
    self._read_reference_information()
    if self.path and blockStart and innerStart:
      self.fh.seek(blockStart,innerStart)
    

  def __iter__(self):
    return self

  def next(self):
    e = self.read_entry()
    if not e:
      raise StopIteration
    else: return e

  def read_entry(self):
    bstart = self.fh.get_block_start()
    innerstart = self.fh.get_inner_start()
    b = self.fh.read(4) # get block size bytes
    if not b: return None
    block_size = struct.unpack('<i',b)[0]
    #print 'block_size '+str(block_size)
    bam = BAM(self.fh.read(block_size),self.ref_names,fileName=self.path,blockStart=bstart,innerStart=innerstart)
    return bam
  def _read_reference_information(self):
    for n in range(self.n_ref):
      l_name = struct.unpack('<i',self.fh.read(4))[0]
      name = self.fh.read(l_name).rstrip('\0')
      l_ref = struct.unpack('<i',self.fh.read(4))[0]
      self.ref_lengths[name] = l_ref
      self.ref_names.append(name)

  def _read_top_header(self):
    magic = self.fh.read(4)
    l_text = struct.unpack('<i',self.fh.read(4))[0]
    self.header_text = self.fh.read(l_text).rstrip('\0')
    self.n_ref = struct.unpack('<i',self.fh.read(4))[0]

def _parse_bam_data_block(bin_in,ref_names):
  v = {}
  data = StringIO(bin_in)
  refID = struct.unpack('<i',data.read(4))[0]
  v['rname'] = ref_names[refID]
  pos = struct.unpack('<i',data.read(4))[0]
  v['pos'] = pos+1
  bin_mq_nl = struct.unpack('<I',data.read(4))[0]
  bin =  bin_mq_nl >> 16 
  mapq = (bin_mq_nl & 0xFF00) >> 8
  v['mapq'] = mapq
  l_read_name = bin_mq_nl & 0xFF
  flag_nc = struct.unpack('<I',data.read(4))[0]
  flag = flag_nc >> 16
  v['flag'] = flag
  n_cigar_op = flag_nc & 0xFF
  l_seq = struct.unpack('<i',data.read(4))[0]
  next_refID = struct.unpack('<i',data.read(4))[0]
  v['rnext'] = ref_names[next_refID]
  next_pos = struct.unpack('<i',data.read(4))[0]
  v['pnext'] = next_pos+1
  tlen = struct.unpack('<i',data.read(4))[0]
  v['tlen'] = tlen
  read_name = data.read(l_read_name).rstrip('\0')
  v['qname'] = read_name
  cigar_bytes = data.read(n_cigar_op*4)
  v['cigar_bytes'] = cigar_bytes
  seq_bytes = data.read((l_seq+1)/2)
  v['seq_bytes'] = seq_bytes
  qual_bytes = data.read(l_seq)
  v['qual_bytes'] = qual_bytes
  extra_bytes = data.read()
  v['extra_bytes'] = extra_bytes
  return v

def _bin_to_qual(qual_bytes):
  return ''.join([chr(struct.unpack('B',x)[0]+33) for x in qual_bytes])
def _bin_to_seq(seq_bytes):
  global _bam_char
  seq = ''.join([''.join([''.join([chr(z+97).translate(_bam_char) for z in  [y>>4,y&0xF]]) for y in struct.unpack('<B',x)]) for x in seq_bytes]).rstrip('=')
  return seq
def _bin_to_cigar(cigar_bytes):
  global _bam_ops
  cigar_packed = [struct.unpack('<I',x)[0] for x in \
             [cigar_bytes[i:i+4] for i in range(0,len(cigar_bytes),4)]]
  cigar_array = [[c >> 4, str(c &0xF).translate(_bam_ops)] for c in cigar_packed]
  cigar_seq = ''.join([''.join([str(x[0]),x[1]]) for x in cigar_array])
  return [cigar_array,cigar_seq]
#Pre all the reamining bytes of an entry
#Post an array of 
# 1. A dict keyed by Tag with {'type':,'value':} where value is a string unless type is i
# 2. A string of the remainder
def _bin_to_extra(extra_bytes):
  global _bam_value_type
  extra = StringIO(extra_bytes)
  tags = {}
  rem = ''
  while extra.tell() < len(extra_bytes):
    tag = extra.read(2)
    val_type = extra.read(1)
    if val_type == 'Z':
      rem += tag+':'
      rem += val_type+':'
      p = re.compile('([!-~])')
      m = p.match(extra.read(1))
      vre = ''
      while m:
        vre += m.group(1)
        c = extra.read(1)
        #print c
        m = p.match(c)
      rem += vre+"\t"
      tags[tag] = {'type':val_type,'value':vre}
    elif val_type == 'A':
      rem += tag+':'
      rem += val_type+':'
      vre = extra.read(1)
      rem += vre+"\t"      
      tags[tag] = {'type':val_type,'value':vre}      
    elif val_type in _bam_value_type:
      rem += tag+':'
      rem += 'i'+':'
      val = struct.unpack(_bam_value_type[val_type][1],extra.read(_bam_value_type[val_type][0]))[0]
      rem += str(val)+"\t"
      tags[tag] = {'type':val_type,'value':val}
    elif val_type == 'B':
      sys.sterr.write("WARNING array not implmented\n")
      continue
      rem += tag+':'
      rem += val_type+':'
      array_type = _bam_value_type[extra.read(1)]
      element_count = struct.unpack('<I',extra.read(4))[0]
      array_bytes = extra.read(element_count*_bam_value_type[array_type][0])
      for by in [array_bytes[i:i+_bam_value_type[array_type][0]] for i in range(0,len(array_bytes),_bam_value_type[array_type][0])]:
        aval = struct.unpack(_bam_value_type[array_type][1],by)
  return [tags,rem.rstrip("\t")]


class BGZF:
  # Methods adapted from biopython's bgzf.py
  def __init__(self,filename,blockStart=None,innerStart=None):
    self.path = filename
    self.fh = open(filename,'rb')
    if blockStart: self.fh.seek(blockStart)
    self._block_start = 0
    #self.pointer = 0
    #holds block_size and data
    self._buffer = self._load_block()
    self._buffer_pos = 0
    if innerStart: self._buffer_pos = innerStart
  def get_block_start(self):
    return self._block_start
  def get_inner_start(self):
    return self._buffer_pos
  def seek(self,blockStart,innerStart):
    self.fh.seek(blockStart)
    self._buffer_pos = 0
    self._buffer = self._load_block()
    self._buffer_pos = innerStart
  def read(self,size):
    done = 0 #number of bytes that have been read so far
    v = ''
    while True:
      if size-done < len(self._buffer['data']) - self._buffer_pos:
        v +=  self._buffer['data'][self._buffer_pos:self._buffer_pos+(size-done)]
        self._buffer_pos += (size-done)
        #self.pointer += size
        return v
      else: # we need more buffer
        vpart = self._buffer['data'][self._buffer_pos:]
        self._buffer = self._load_block()
        v += vpart
        self._buffer_pos = 0
        if len(self._buffer['data'])==0: return v
        done += len(vpart)

  def _load_block(self):
    #pointer_start = self.fh.tell()
    if not self.fh: return {'block_size':0,'data':''}
    self._block_start = self.fh.tell()
    magic = self.fh.read(4)
    if len(magic) < 4:
      #print 'end?'
      #print len(self.fh.read())
      return {'block_size':0,'data':''}
    gzip_mod_time, gzip_extra_flags, gzip_os,extra_len = struct.unpack("<LBBH",self.fh.read(8))
    pos = 0
    block_size = None
    #get block_size
    while pos < extra_len:
      subfield_id = self.fh.read(2)
      subfield_len = struct.unpack("<H",self.fh.read(2))[0]
      subfield_data = self.fh.read(subfield_len)
      pos += subfield_len+4
      if subfield_id == 'BC':
        block_size = struct.unpack("<H",subfield_data)[0]+1
    #block_size is determined
    deflate_size = block_size - 1 - extra_len - 19
    d = zlib.decompressobj(-15)
    data = d.decompress(self.fh.read(deflate_size))+d.flush()
    expected_crc = self.fh.read(4)
    expected_size = struct.unpack("<I",self.fh.read(4))[0]
    if expected_size != len(data):
      sys.stderr.write("ERROR unexpected size\n")
      sys.exit()
    crc = zlib.crc32(data)
    if crc < 0:  crc = struct.pack("<i",crc)
    else:  crc = struct.pack("<I",crc)
    if crc != expected_crc:
      sys.stderr.write("ERROR crc fail\n")
      sys.exit()
    return {'block_size':block_size, 'data':data}
    #print crc
    #print len(data)
    #print expected_crc
    #print expected_size
    #print 'hello'

#Some helper functions for reading binary
