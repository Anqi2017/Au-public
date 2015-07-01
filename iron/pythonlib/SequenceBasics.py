import re, sys, hashlib
from FileBasics import GenericFileReader


class GenericFastqFileReader:
  def __init__(self,filename):
    self.filename = filename
    self.gfr = GenericFileReader(self.filename)
    self.previous_name = None

  def close(self):
    self.gfr.close()

  def read_entry(self):
    line1 = self.gfr.readline()
    if not line1:
      return False
    line2 = self.gfr.readline()
    if not line2:
      sys.stderr.write("Warning: Improperly terminated fastq file line count not a multiple of 4\n")
    line3 = self.gfr.readline()
    if not line3:
      sys.stderr.write("Warning: Improperly terminated fastq file line count not a multiple of 4\n")
    line4 = self.gfr.readline()
    if not line4:
      sys.stderr.write("Warning: Improperly terminated fastq file line count not a multiple of 4\n")
    m = re.match('^@([^\t]+)',line1.rstrip())
    if not m:
      sys.stderr.write("Warning: Could not read name\n")
    t = {}
    t['name'] = m.group(1)
    t['seq'] = line2.rstrip()
    t['quality'] = line4.rstrip()
    return t

class GenericFastaFileReader:
  def __init__(self,filename):
    self.filename = filename
    self.gfr = GenericFileReader(self.filename)
    self.previous_name = None
  def close(self):
    self.gfr.close()
  def read_entry(self):
    buffer = ''
    original = ''
    t = {}
    t['name'] = ''
    t['seq'] = ''
    t['original'] = ''
    while True:
      newline = self.gfr.readline()
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

# an upgrade to the old sequence_basics set

    
# pre: A nucleotide sequence.  no spaces no whitespace
# post Reverse complemented sequence
def rc(seq):
  if re.search('[uU]',seq) and re.search('[tT]',seq):
    print "Mix of Uu and Tt in sequence.  I don't know what it is."
    sys.exit()
  rna = 0
  if re.search('[uU]',seq):
    rna = 1
  o = ''
  for i in range(0,len(seq)):
    c = seq[i]
    if c == 'A' and rna == 0: o = 'T' + o
    elif c == 'A' and rna == 1: o = 'U' + o
    elif c == 'C': o = 'G' + o
    elif c == 'G': o = 'C' + o
    elif c == 'T': o = 'A' + o
    elif c == 'U': o = 'A' + o
    elif c == 'a' and rna == 0: o = 't' + o
    elif c == 'a' and rna == 1: o = 'u' + o
    elif c == 'c': o = 'g' + o
    elif c == 'g': o = 'c' + o
    elif c == 't': o = 'a' + o
    elif c == 'u': o = 'a' + o
    else: o = c + o  #just use this odd character as is.
  return o


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
