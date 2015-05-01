import GenePredBasics, SequenceBasics, PSLBasics
from FileBasics import GenericFileReader
import re, sys, os
import subprocess

class PSLtoSAMconversionFactory:
  # Based on the 3 Mar 2015 Sam Specification
  # Can take a lot of RAM because of needing to store the fasta
  def __init__(self):
    self.reads = {}
    self.qualities = {}
    self.min_intron_size = 68
    self.reads_set = False
    self.qualities_set = False
    self.mapping_counts_set = False
    self.ref_genome_set = False
    self.skip_directionless_splice = False
    self.set_canon()
    self.set_revcanon()
  def set_skip_directionless_splice(self):
    self.skip_directionless_splice = True
  def set_canon(self):
    v = set()
    v.add('GT-AG')
    v.add('GC-AG')
    v.add('AT-AC')
    self.canonical = v
  def set_revcanon(self):
    v = set()
    v.add('CT-AC')
    v.add('CT-GC')
    v.add('GT-AT')
    self.revcanonical = v
  def set_mapping_counts(self,psl_filename):
    self.mapping_counts_set = True
    gfr0 = GenericFileReader(psl_filename)
    qcnts = {}
    while True:
      line = gfr0.readline()
      if not line: break
      psle = PSLBasics.line_to_entry(line.rstrip())
      if psle['qName'] not in qcnts: qcnts[psle['qName']] = 0
      qcnts[psle['qName']] += 1
    gfr0.close()
    self.mapping_counts = qcnts

  def set_min_intron_size(self,intron_size):
    self.min_intron_size = intron_size

  def set_reference_genome(self,ref_genome):
    self.ref_genome_set = True
    self.ref_genome = SequenceBasics.read_fasta_into_hash(ref_genome)

  def convert_line(self,psl_line):
    pe = PSLBasics.line_to_entry(psl_line)
    #work on the positive strand case first
    cigar = '*'
    blocks = len(pe['blockSizes'])
    starts = pe['qStarts']
    #if pe['strand'] == '-':
    #  starts = [x for x in reversed(pe['qStarts_actual'])]
    #  print 'isrev'
    q_coord_start = starts[0]+1 # base-1 converted starting position
    q_coord_end = starts[blocks-1]+pe['blockSizes'][blocks-1] # base-1 position
    t_coord_start = pe['tStarts'][0]+1 # base-1 converted starting position
    t_coord_end = pe['tStarts'][blocks-1]+pe['blockSizes'][blocks-1] # base-1 position
    if pe['qName'] not in self.reads and self.reads_set is True:
      sys.stderr.write("Warning: qName "+pe['qName']+" was not found in reads\n")
    # we will clip the query sequence to begin and end from the aligned region
    #q_seq = ''
    #if self.reads_set:
    #  q_seq = self.reads[pe['qName']]

    # 1. Get the new query to output
    q_seq_trimmed = '*'
    if self.reads_set:
      q_seq_trimmed = self.reads[pe['qName']]
      if pe['strand'] == '-':
        q_seq_trimmed = SequenceBasics.rc(q_seq_trimmed)
      q_seq_trimmed = q_seq_trimmed[q_coord_start-1:q_coord_end]

    qual_trimmed = '*'
    if self.qualities_set:
      qual_trimmed = self.qualities[pe['qName']]
      if pe['strand'] == '-':
        qual_trimmed = qual_trimmed[::-1]
      qual_trimmed = qual_trimmed[q_coord_start-1:q_coord_end]
    # 2. Get the cigar string to output
    prev_diff = t_coord_start-q_coord_start
    cigar = ''
    #for i in range(0,blocks):
    #  current_diff = pe['tStarts'][i]-starts[i]
    #  delta = current_diff - prev_diff
    #  #print delta
    #  if delta >= self.min_intron_size:
    #    cigar += str(abs(delta))+'N'
    #  elif delta > 0: # we have a
    #    cigar += str(abs(delta))+'D'
    #  elif delta < 0: # we have a
    #    cigar += str(abs(delta))+'I'
    #  cigar += str(pe['blockSizes'][i])+'M' # our matches
    #  #print current_diff
    #  prev_diff = current_diff
    qstarts = [x-pe['qStarts'][0] for x in pe['qStarts']]
    tstarts = [x-pe['tStarts'][0] for x in pe['tStarts']]
    query_index = 0
    target_index = 0
    junctions = []
    for i in range(0,blocks):
      qdif = qstarts[i] - query_index
      tdif = tstarts[i] - target_index
      if qdif > 0:  # we have to insert
        cigar += str(qdif) + 'I'
      if tdif > self.min_intron_size: # we have an intron
        cigar += str(tdif) + 'N'
        junctions.append(i)
      elif tdif > 0: # we have to delete
        cigar += str(tdif) + 'D'
      cigar += str(pe['blockSizes'][i]) + 'M'
      query_index = qstarts[i]+pe['blockSizes'][i]
      target_index = tstarts[i]+pe['blockSizes'][i]
    ### cigar done
    # inspect junctions if we have a ref_genome
    spliceflag_set = False
    if self.ref_genome_set:
      canon = 0
      revcanon = 0
      for i in junctions: #blocks following a junction
        left_num = pe['tStarts'][i-1]+pe['blockSizes'][i-1]
        left_val = self.ref_genome[pe['tName']][left_num:left_num+2].upper()
        right_num = pe['tStarts'][i-1]-2
        right_val = self.ref_genome[pe['tName']][right_num:right_num+2].upper()
        junc = left_val + '-' + right_val
        if junc in self.canonical: canon += 1
        if junc in self.revcanonical: revcanon += 1
      if canon > revcanon: 
        spliceflag_set = True
        spliceflag = '+'
      elif revcanon > canon:
        spliceflag_set = True
        spliceflag = '-'
    # if we have junctions, and we should be setting direction but 
    # we can't figure out the direction skip ambiguous direction
    if len(junctions) > 0 and self.skip_directionless_splice and spliceflag_set == False:
      return False
    samline =  pe['qName'] + "\t"        # 1. QNAME
    if pe['strand'] == '-':
      samline += '16' + "\t"             # 2. FLAG
    else:
      samline += '0' + "\t"
    samline += pe['tName'] + "\t"        # 3. RNAME
    samline += str(t_coord_start) + "\t" # 4. POS
    samline += '0' + "\t"                # 5. MAPQ
    samline += cigar + "\t"         # 6. CIGAR
    samline += '*' + "\t"           # 7. RNEXT
    samline += '0' + "\t"           # 8. PNEXT
    samline += '0' + "\t"           # 9. TLEN
    samline += q_seq_trimmed + "\t" # 10. SEQ
    samline += qual_trimmed + "\t"  # 11. QUAL
    if spliceflag_set:
      samline += 'XS:A:'+spliceflag + "\t"
    if self.ref_genome_set:
      samline += 'NH:i:'+str(self.mapping_counts[pe['qName']]) + "\t"
    samline += 'XC:i:'+str(len(junctions)) + "\t"
    samline += 'NM:i:0'
    return samline

  def set_read_fasta(self,read_fasta_file):
    self.reads_set = True
    gfr = SequenceBasics.GenericFastaFileReader(read_fasta_file)
    self.reads = {}
    while True:
      e = gfr.read_entry()
      if not e: break
      if e['name'] in self.reads:
        sys.stderr.write("Warning duplicate name in fasta file, could be big problems on sequence assignment.\n")
      self.reads[e['name']] = e['seq'].upper()
    gfr.close()
    return
  def set_read_fastq(self,read_fastq_file):
    self.reads_set = True
    self.qualities_set = True
    gfr = SequenceBasics.GenericFastqFileReader(read_fastq_file)
    self.reads = {}
    self.qualities = {}
    while True:
      e = gfr.read_entry()
      if not e: break
      if e['name'] in self.reads:
        sys.stderr.write("Warning duplicate name in fasta file, could be big problems on sequence assignment.\n")
      self.reads[e['name']] = e['seq'].upper()
      self.qualities[e['name']] = e['quality']
    gfr.close()
    return
     

def construct_header_from_reference_fasta(ref_fasta_filename):
  g = SequenceBasics.read_fasta_into_hash(ref_fasta_filename)
  chrs = {}
  for name in sorted(g):
    chrs[name] = len(g[name])
    sys.stderr.write(name+" is there at length "+str(len(g[name]))+"\n")
  header = ''
  header += "@HD\tVN:1.0\tSO:coordinate\n"
  for chr in sorted(chrs):
    header += "@SQ\tSN:"+chr+"\tLN:"+str(chrs[chr])+"\n"
  header += "@PG\tID:SamBasics.py\tVN:1.0\n"
  return header 



# Pre: Requires an indexed bam file
# 
class RandomAccessCoordinateReader:
  def __init__(self,bam_filename,chromosome,start,finish):
    self.filename = bam_filename
    if not re.search('\.bam$',self.filename): 
      sys.stderr.write("Error: not bam file extension.\n")
      sys.exit()
    if not os.path.isfile(self.filename+'.bai'):
      sys.stderr.write("Error: no index bam.bai file present.\n")
      sys.exit()
    cmd = 'samtools view '+self.filename+' '+chromosome+':'+str(start)+'-'+str(finish)
    args = cmd.split()
    self.process = subprocess.Popen(args,bufsize=0,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  def close(self):
    if self.process:
      self.process.kill()

  def readline(self):
      return self.process.stdout.readline()

  def readentry(self):
    line = self.process.stdout.readline()
    if not line: return line
    return sam_line_to_dictionary(line.rstrip())

# Pre: A sam entry (dictionary)
# Post: a bed file line describing the start and stop
def entry_to_blocked_bed(entry,color):
  if entry['rname'] == '*':  return False
  ostring = entry['rname'] + "\t" + str(entry['pos']-1) + "\t"
  #print entry
  #print ostring
  z = entry['pos']-1
  block_count = 0
  block_starts = []
  block_sizes = []
  for c in entry['cigar_array']:
    # entry maps
    if re.match('[MISX=]',c['op']):
      # here is where we should output
      block_starts.append(z-(entry['pos']-1))
      block_sizes.append(c['val'])
      z += c['val']
      block_count += 1
    # entry is a gap
    if re.match('[DNH]',c['op']):
      z+= c['val']
  endfeature = block_starts[block_count-1]+block_sizes[block_count-1]+(entry['pos']-1)
  ostring += str(endfeature) + "\t" # chromEnd
  ostring += entry['qname'] + "\t" # name
  ostring += "1" + "\t" # score
  strand = get_entry_strand(entry)
  ostring += strand + "\t" # strand
  ostring += str(entry['pos']-1) + "\t" #thickStart
  ostring += str(endfeature) + "\t"
  ostring += color + "\t" #itemRgb
  ostring += str(block_count) + "\t" # block count
  ostring += ",".join([str(x) for x in block_sizes])+"," + "\t"   #blockSizes
  ostring += ",".join([str(x) for x in block_starts])+","  #blockStarts
  #ostring += 
  return ostring

def get_entry_strand(entry):
  if check_flag(entry['flag'],16):
    return '+'
  else:
    return '-'
  #print entry['remainder']
  #if re.search('XS:A:+',entry['remainder']): return '+'
  #elif re.search('XS:A:-',entry['remainder']): return '-'
  #else: 
  #  sys.stderr.write("Error did not find strand information for "+entry['rname']+"\n")
  #  sys.exit()
  #return False

class GenericSamReader:
  def __init__(self,filename):
    self.filename = filename
    if re.search('\.bam$',self.filename): # do it as a gzipped stream
      cmd = 'samtools view '+filename
      args = cmd.split()
      self.process = subprocess.Popen(args,stdout=subprocess.PIPE)
    else:
      cmd = 'samtools view -S '+filename
      args = cmd.split()
      self.process = subprocess.Popen(args,bufsize=0,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

  def close(self):
    if self.process:
      self.process.kill()

  def readline(self):
      return self.process.stdout.readline()

# pre: sam file name, genepred name 
# post: list of coordinates in the reference format
#       these coordinates are zero indexed for both the start and end coordinate
#       note thats different than gpd or psl or even wig or bed
#       <read name> <genepred entry name> <chromosome:coord1-coord2,coord3-coord4,...>
def convert_directionless_gpd_alignment_to_reference(sam_filename,genepred_filename,out_map):
  conv = GenePredBasics.get_directionless_gpd_conversion(genepred_filename)
  ofile = open(out_map,'w')
  with open(sam_filename) as samfile:
    for line in samfile:
      line = line.rstrip()
      if re.match('^@[A-Z][A-Z]\s',line): continue #skip header
      d = sam_line_to_dictionary(line)
      if d['rname'] == '*': continue #skip unmapped
      startposition = d['pos']-1
      readcoord = []
      z = 0
      for entry in d['cigar_array']:
        if re.match('[MISX=]',entry['op']):  # all the entries that map to the read
          for i in range(0,entry['val']):
            if re.match('[M=X]',entry['op']): #all the entries that match the reference alignment
              readcoord.append(conv[d['rname']]['coordinates'][startposition+z])
              z+=1
            # lets ignore insertions for now
            #else:
            #  readcoord.append('*')
        if re.match('[DNH]',entry['op']):
          z+= entry['val']      
      abbrev = conv[d['rname']]['chrom']+':'+SequenceBasics.collapse_coordinate_array(readcoord)
      ofile.write(d['qname'] + "\t" + d['rname'] + "\t" + abbrev + "\n")
  ofile.close()

# pre:       A line from a sam file
# post:      a string with the coordiantes of the alignment

def get_coordinates(sam_line):
  f = sam_line.rstrip().split("\t")
  name = f[0]
  coordinate = ''
  if f[2] == '*':
    return [name, coordinate]
  coordinate = f[2]+':'+str(f[3])+':'+f[5]
  return [name,coordinate]


#pre: a flag from a sam file, in integer format
#     a bit to convert, given as a hex number ie 0x10
#post: returns true if the flag is set on
def is_header(line):
  if re.match('^@',line):
    f = line.rstrip().split("\t")
    if(len(f) > 9):
      return False
    return True
  return False


#pre: a flag from a sam file, in integer format
#     a bit to convert, given as a hex number ie 0x10
#post: returns true if the flag is set on
def check_flag(flag,inbit):
  if flag & inbit:
    return True
  return False

#pre: a line from a sam file that is not a header entry
#post: a dictionary with entries named like the manual
def sam_line_to_dictionary(line):
  f = line.rstrip().split("\t")
  d = {}
  d['qname'] = f[0]
  d['flag'] = int(f[1])
  d['rname'] = f[2]
  d['pos'] = int(f[3])
  d['mapq'] = f[4]
  cigar = parse_cigar(f[5])
  d['cigar_array'] = cigar
  d['cigar'] = f[5]
  d['rnext'] = f[6]
  d['pnext'] = int(f[7])
  d['tlen'] = int(f[8])
  d['seq'] = f[9]
  d['qual'] = f[10]
  d['remainder'] = ''
  if len(f) > 11:
    for i in range(11,len(f)):
      d['remainder'] += f[i]+" "
    d['remainder'] = d['remainder'].rstrip(" ")
  return d

# pre: CIGAR string
# post: an array of cigar string entries
def parse_cigar(cigar):
  v = re.findall('([0-9]+[A-Z])',cigar)
  vals = []
  for val in v:
    m = re.match('(\d+)([A-Z])',val)
    d = {}
    d['op'] = m.group(2)
    d['val'] = int(m.group(1))
    vals.append(d)
  return vals

# index 1 coordinates
def get_base_at_coordinate(entry,chr,coord):
  if entry['rname'] != chr:
    return False
  #print chr + "\t" + str(coord)
  z = entry['pos']
  bases = list(entry['seq'])
  b = 0
  for c in entry['cigar_array']:
    # entry maps
    if re.match('[MISX=]',c['op']):
      # here is where we should output
      for i in range(0,c['val']):
        if int(z) == int(coord):
          return bases[b]
        b += 1
        z += 1
    # entry is a gap
    if re.match('[DNH]',c['op']):
      z+= c['val']
  return False
