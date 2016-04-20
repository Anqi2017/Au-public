import struct, sys, zlib, StringIO, time

class reader:
  # Methods adapted from biopython's bgzf.py
  def __init__(self,handle,blockStart=None,innerStart=None):
    self.fh = handle
    self._pointer = 0
    self._block_start = 0
    if blockStart: 
      self.fh.seek(blockStart)
      self._pointer = blockStart
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
    self._pointer = blockStart
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
    #self._block_start = self.fh.tell()
    self._block_start = self._pointer
    magic = self.fh.read(4)
    self._pointer += 4
    if len(magic) < 4:
      #print 'end?'
      #print len(self.fh.read())
      return {'block_size':0,'data':''}
    gzip_mod_time, gzip_extra_flags, gzip_os,extra_len = struct.unpack("<LBBH",self.fh.read(8))
    self._pointer += 8
    pos = 0
    block_size = None
    #get block_size
    while pos < extra_len:
      subfield_id = self.fh.read(2)
      self._pointer += 2
      subfield_len = struct.unpack("<H",self.fh.read(2))[0]
      self._pointer += 2
      subfield_data = self.fh.read(subfield_len)
      pos += subfield_len+4
      if subfield_id == 'BC':
        block_size = struct.unpack("<H",subfield_data)[0]+1
        #print 'blocksize :'+str(block_size)
    #block_size is determined
    deflate_size = block_size - 1 - extra_len - 19
    #deflate_size = block_size - extra_len - 19
    d = zlib.decompressobj(-15)
    data = d.decompress(self.fh.read(deflate_size))+d.flush()
    self._pointer += deflate_size
    expected_crc = self.fh.read(4)
    self._pointer += 4
    expected_size = struct.unpack("<I",self.fh.read(4))[0]
    self._pointer += 4
    #print len(data)
    #print expected_size
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

class writer:
  #  Give it the handle of the stream to write to
  def __init__(self,handle):
    #self.path = filename
    self.fh = handle
    self.buffer_size = 64000
    self.buffer = bytearray()
  def __del__(self): 
    self.close()
  def write(self,bytes):
    self.buffer+=bytes
    if len(self.buffer) < self.buffer_size:
      return True
    while len(self.buffer) >= self.buffer_size:
      dobytes = self.buffer[0:self.buffer_size]
      self.buffer = self.buffer[self.buffer_size:]
      self._do_block(dobytes)
    return

  def close(self):
    if len(self.buffer) == 0:  return True
    self._do_block(self.buffer)
    self.buffer = []
    return True

  def _do_block(self,bytes):
    # now we can output this  
    isize = len(bytes)
    s = StringIO.StringIO(bytes)
    d = zlib.compressobj(9,zlib.DEFLATED,-zlib.MAX_WBITS)
    data = d.compress(str(bytes))+d.flush()
    datasize = len(data)
    output = bytearray()
    output += struct.pack('<B',31)  #IDentifier1
    output += struct.pack('<B',139) #IDentifier2
    output += struct.pack('<B',8)   #Compression Method
    output += struct.pack('<B',4)   #FLaGs
    output += struct.pack('<I',int(time.time())) #Modification TIME
    output += struct.pack('<B',0)   #eXtra FLags
    output += struct.pack('<B',0x03)   #Operating System = Unix
    output += struct.pack('<H',6)   #eXtra LENgth
    # Subfields
    output += struct.pack('<B',66) #Subfield Identifier 1
    output += struct.pack('<B',67) # Subfield Identifier 2
    output += struct.pack('<H',2) #Subfield Length
    outsize = datasize+19+6
    output += struct.pack('<H',outsize) #Total block size minus one
    output += data
    crc = zlib.crc32(str(bytes))
    if crc < 0:  output += struct.pack("<i",crc)
    else:  output+= struct.pack("<I",crc)
    output += struct.pack("<I",isize) #isize
    self.fh.write(str(output))    
