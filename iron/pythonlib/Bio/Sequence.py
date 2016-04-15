import re, sys, hashlib, base64, gzip, os
from string import maketrans
#from multiprocessing.sharedctypes import RawArray

class FastaHandleReader:
  def __init__(self,in_handle):
    self.fh = in_handle
    self.previous_name = None
  def close(self):
    self.fh.close()
  def __iter__(self):
    return self
  def next(self):
    v = self.read_entry()
    if not v:
      raise StopIteration
    else:
      return v
  def read_entry(self):
    buffer = ''
    original = ''
    t = {}
    t['name'] = ''
    t['seq'] = ''
    t['original'] = ''
    while True:
      newline = self.fh.readline()
      if not self.previous_name and not newline:
        # no name in the buffer and new data being input, exit
        return None
      if not newline:
        # end of the line, then finish it
        t['name'] = self.previous_name
        t['seq'] = buffer
        t['original'] = original
        self.previous_name = None
        t['original'] = '>'+t['name'] + "\n" + t['original']
        return t
      m = re.match('^>(.*)$',newline.rstrip())
      if not self.previous_name and m:
        self.previous_name = m.group(1)
        #special case of our first entry
        continue
      if m:
        t['name'] = self.previous_name
        t['seq'] = buffer
        t['original'] = original
        self.previous_name = m.group(1)
        t['original'] = '>'+t['name'] + "\n" + t['original']
        return t
      buffer += newline.rstrip()
      original += newline

class FastqHandleReader:
  def __init__(self,in_handle):
    self.fh = in_handle
  def close(self):
    self.fh.close()
  def read_entry(self):
    t = {}
    t['name'] = ''
    t['seq'] = ''
    t['original'] = ''
    t['qual'] = ''
    line1 = self.fh.readline()
    if not line1: return None
    line2 = self.fh.readline()
    if not line2: return None
    line3 = self.fh.readline()
    if not line3: return None
    line4 = self.fh.readline()
    if not line4: return None
    # end of the line, then finish it
    m1 = re.match('^@(.*)$',line1.rstrip())
    if not m1:
      sys.stderr.write('ERROR: '+str(line1)+"\n")
      sys.exit()
    t['name'] = m1.group(1)
    t['seq'] = line2.rstrip()
    t['qual'] = line4.rstrip()
    t['original'] = line1+line2+line3+line4
    return t
# an upgrade to the old sequence_basics set

    
def rc(seq):
  complement = maketrans('ACTGUNXactgunx','TGACANXtgacanx')
  return seq.translate(complement)[::-1]
# pre:  A fasta file name (input)
# post: Dictionary containing names and counts
#       names only the first non-whitespace, not the whole header
# Modifies: file IO

def counts_by_name_from_fasta(in_fasta):
  reads = {}
  with open(in_fasta) as f:
    for line in f:
      line = line.rstrip()
      m = re.match('^>([\S]+)',line)
      if m:
        if not m.group(1) in reads:
          reads[m.group(1)] = 0
        reads[m.group(1)]+=1
  return reads

# pre:  A fasta file name (input), a fasta file name (output)
# post: Writes to the fasta file name (output)
#       Header information is destroyed and a new unqiue indicator in its place
# Modifies: file IO

def fasta_to_unique_name_fasta(in_fasta,out_fasta):
  ofile=open(out_fasta,'w')
  i = 0
  with open(in_fasta) as f:
    for line in f:
      line = line.rstrip()
      if re.match('^>',line):
        i+=1
        ofile.write('>SR'+str(i)+"\n")
      else:
        ofile.write(line+"\n")
  ofile.close()

# pre:  A fastq file name (input), a fastq file name (output)
# post: Writes to the fastq file name (output)
#       Header information is lost and a new unique indicator in its place
# Modifies: file IO
def fastq_to_unique_name_fastq(in_fastq,out_fastq):
  ofile=open(out_fastq,'w')
  i = 0
  with open(in_fastq) as f:
    while True:
      l1 = f.readline().rstrip()
      if not l1: 
        break
      l2 = f.readline().rstrip()
      l3 = f.readline().rstrip()
      l4 = f.readline().rstrip()
      i+=1
      ofile.write('@SR'+str(i)+"\n")
      ofile.write(l2+"\n")
      ofile.write('+SR'+str(i)+"\n")
      ofile.write(l4+"\n")
  ofile.close()

# pre:  A fastq file name (input), a list of read names, a fastq file name (output)
# post: Writes to the fastq file name (output)
#       only writes ones where the read name matches, matches on first nonwhitespace
# Modifies: file IO
def write_fastq_subset(in_fastq,names,out_fastq):
  nameset = set()
  for name in names: nameset.add(name)
  ofile=open(out_fastq,'w')
  with open(in_fastq) as f:
    while True:
      l1 = f.readline().rstrip()
      if not l1: 
        break
      l2 = f.readline().rstrip()
      l3 = f.readline().rstrip()
      l4 = f.readline().rstrip()
      m = re.match('^@([\S]+)')
      name = ''
      if m: name = m.group(1)
      if name in nameset:
        ofile.write(l1+"\n")
        ofile.write(l2+"\n")
        ofile.write(l3+"\n")
        ofile.write(l4+"\n")
  ofile.close()

# pre:  A fastq file name (input)
# post: Dictionary containing read names and the number of times present
#       name is the first non-whitespace, not the whole header
# Modifies: none
def counts_by_name_from_fastq(in_fastq,out_fastq):
  reads = {}
  with open(in_fastq) as f:
    while True:
      l1 = f.readline().rstrip()
      if not l1: 
        break
      l2 = f.readline().rstrip()
      l3 = f.readline().rstrip()
      l4 = f.readline().rstrip()
      m = re.match('@([\S]+)')
      if m:
        if not m.group(1) in reads:
          reads[m.group(1)] = 0
        reads[m.group(1)] += 1
  return reads

# pre:       A fasta file name
# post:      An array of dictionaries,
#            with 'name' and 'seq' entires.
# modifies:  none

def write_fasta_subset(in_fasta,names,out_fasta):
  nameset = set()
  for name in names: nameset.add(name)
  ofile = open(out_fasta,'w')
  with open(in_fasta) as f:
    seq = ''
    head = ''
    prog = re.compile('^>(.+)')
    for line in f:
      line = line.rstrip()
      m = prog.match(line)
      if m:
        if seq != '':
          mhead = re.match('^(\S+)',head)
          usename = ''
          if mhead: usename = mhead.group(1)
          if usename in nameset:
            ofile.write(">"+str(head)+"\n"+seq+"\n")
          seq = ''
        head = m.group(1)
      else:
         seq = seq + line
    if seq != '':
      mhead = re.match('^(\S+)',head)
      usename = ''
      if mhead: usename = mhead.group(1)
      if usename in nameset:
        ofile.write(">"+str(head)+"\n"+seq+"\n")

# pre:       A fasta file name
# post:      An array of dictionaries,
#            with 'name' and 'seq' entires.
# modifies:  none

def read_fasta_into_array(fasta_filename):
  seqs = []
  with open(fasta_filename) as f:
    seq = ''
    head = ''
    prog = re.compile('^>(.+)')
    for line in f:
      line = line.rstrip()
      m = prog.match(line)
      if m:
        if seq != '':
          val = {}
          val['name'] = head
          val['seq'] = seq
          seqs.append(val)
          seq = ''
        head = m.group(1)
      else:
         seq = seq + line
    if seq != '':
      val = {}
      val['name'] = head
      val['seq'] = seq
      seqs.append(val)
  return seqs

# pre:       A fasta file name
# post:      A dictionary of sequences keyed by name
# modifies:  none

def read_fasta_into_hash(fasta_filename):
  seqs = {}
  with open(fasta_filename) as f:
    seq = ''
    head = ''
    prog = re.compile('^>(.+)')
    for line in f:
      line = line.rstrip()
      m = prog.match(line)
      if m:
        if seq != '':
          seqs[head] = seq
          seq = ''
        head = m.group(1)
      else:
         seq = seq + line
    if seq != '':
      seqs[head] = seq
  return seqs

# pre: A coordiante array
# post: a more bed format like string
def collapse_coordinate_array(readcoords):
  sc = readcoords[0]
  lc = readcoords[0] 

  # readcoodinates should be size of read length unless a sequnece in the read didn't map (was an insertion)
  cstring = ''
  for i in range(1,len(readcoords)):
    if readcoords[i] > int(lc)+1:
      cstring = cstring + ',' + str(sc) + "-" + str(lc)
      sc = readcoords[i]
    lc = readcoords[i]
  cstring = cstring + ',' +  str(sc) + "-" + str(lc)
  cstring = cstring.lstrip(',')
  return cstring

def encode_name(conversion_string):
  compressed_string = zlib.compress(conversion_string,9)
  enc_string = base64.b32encode(compressed_string)
  return 'SZ_'+enc_string.rstrip('=')

def decode_name(safename):
  frag = safename.lstrip('SZ_')
  padding = (8-(len(frag) % 8)) % 8
  c = base64.b32decode(frag+'='*padding)
  return  zlib.decompress(c)

# Slicable fast fasta
class FastaData:
  def __init__(self,data=None,file=None,dict=None):
    self._lengths = {}
    self._seqs = {}
    self._names = []
    # now populate
    if data:
      self._scan_data(data)
    elif file and re.search('\.gz$',file):
      self._scan_data(gzip.open(file,'rb').read())
    elif file:
      self._scan_data(open(file,'rb').read())
    elif dict:
      for name in dict.get_names(): 
        self._names.append(name)
        self._seqs = dict
        self._lengths[name] = len(dict[name])

  def __getitem__(self,key):
    return self._seqs[key]

  def get_sequence(self,chr=None,start=None,end=None,dir=None,rng=None):
    if rng: 
      chr = rng.chr
      start = rng.start
      end = rng.end
      dir = rng.direction
    if not start: start = 1
    if not end: end = self.fai[chr]['length']
    if not dir: dir = '+'
    if dir == '-':
      return rc(self._seqs[chr][start-1:end])
    return self._seqs[chr][start-1:end]

  def _scan_data(self,dat):
    p = re.compile('>([^\n]+)\n([^>]+)')
    pos = 0
    for m in p.finditer(dat):
      self._names.append(m.group(1))
      seq = m.group(2).replace("\n",'')
      self._lengths[m.group(1)] = len(seq)
      self._seqs[m.group(1)] = seq

# Do random access with an indexed Fasta File
# Creates the index if its not there already
# Pre: An uncompressed fasta file
#      Can be called by chromosome and location slices
#          Slices are same as array - zero indexed
# Post: Makes index if doesn't exist upon being called.
#       Can access sequence
# Modifies: File IO reads the fasta, and writes a fasta index file
class FastaFile:
  def __init__(self,fname,index=None):
    self.fname = fname
    self.index = index
    self.fai = {}
    if not self.index:
      if os.path.isfile(self.fname+'.fai'): self.index = self.fname+'.fai'
    if not self.index:
      sys.stderr.write("Warning no index trying to create\n")
      self._make_index()
    self._read_index()
    self.fh = open(fname)

  def __getitem__(self,key):
    chr = FastaFile.Chromosome(self,key)
    if chr.sliced: return chr
    return self.get_sequence(key)

  class Chromosome:
    def __init__(self,outer,chr):
      self.outer = outer
      self.chr = chr
      self.sliced = False

    def __getitem__(self,val):
      self.sliced = True
      print val
      if val.step:  
        sys.stderr.write("ERROR: FastaFile doesn't support step access\n")
        sys.exit()
      clen = self.outer.fai[self.chr]['length']
      return self.outer.get_sequence(self.chr,val.start+1,min(val.stop,clen))

    def __len__(self):
      return self.outer.fai[self.chr]['length']

  def get_sequence(self,chr=None,start=None,end=None,dir=None,rng=None):
    if rng: 
      chr = rng.chr
      start = rng.start
      end = rng.end
      dir = rng.direction
    if not start: start = 1
    if not end: end = self.fai[chr]['length']
    if not dir: dir = '+'
    [sblocks,srem] = divmod(start,self.fai[chr]['linebases'])
    missing_start = sblocks*(self.fai[chr]['linewidth']-self.fai[chr]['linebases'])
    [eblocks,erem] = divmod(end,self.fai[chr]['linebases'])
    missing_end = eblocks*(self.fai[chr]['linewidth']-self.fai[chr]['linebases'])
    pos_start = self.fai[chr]['offset']+start+missing_start-1
    pos_end = self.fai[chr]['offset']+end+missing_end
    self.fh.seek(pos_start)
    v = self.fh.read(pos_end-pos_start).replace("\n",'')
    if dir == '-':
      return rc(v)
    return v

  def _read_index(self):
    with open(self.index) as inf:
      for line in inf:
        v = line.rstrip().split("\t")
        self.fai[v[0]] = {}
        self.fai[v[0]]['name'] = v[0]
        self.fai[v[0]]['length'] = int(v[1])
        self.fai[v[0]]['offset'] = int(v[2])
        self.fai[v[0]]['linebases'] = int(v[3])
        self.fai[v[0]]['linewidth'] = int(v[4])

  def _make_index(self):
    of = None
    try:
      of = open(self.fname+'.fai','w')
    except IOError:
      sys.stderr.write("ERROR: could not open file\n")
      sys.exit()
    self.index = self.fname+'.fai'
    # now we can parse the file
    p1 = re.compile('>([^\r\n]+)([\n\r]+)([^>]+)')
    p2 = re.compile('([^\r\n]+)($|[\n\r]+)')
    pos = 0
    for m1 in p1.finditer(open(self.fname).read()):
      name = m1.group(1)
      pos += len(m1.group(1))+1+len(m1.group(2))
      linewidth_bases = None
      linewidth_bytes = None
      seqlen = 0
      nextoffset = len(m1.group(3))
      for m2 in p2.finditer(m1.group(3)):
        if not linewidth_bases: linewidth_bases = len(m2.group(1))
        elif linewidth_bases != len(m2.group(1)) and len(m2.group(2)) > 0:
          sys.stderr.write("ERROR: irregular line breaks\n")
          sys.exit()
        if not linewidth_bytes: linewidth_bytes = len(m2.group(1))+len(m2.group(2))
        elif linewidth_bytes != len(m2.group(1))+len(m2.group(2)) and len(m2.group(2)) > 0:
          sys.stderr.write("ERROR: irregular bytes line breaks\n")
          sys.exit() 
        seqlen += len(m2.group(1))      
      of.write(name + "\t" + str(seqlen) + "\t"+str(pos)+"\t" + str(linewidth_bases) + "\t" + str(linewidth_bytes)+"\n")
      pos += nextoffset
    of.close()
