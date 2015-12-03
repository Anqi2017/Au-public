import re, sys
import subprocess
import FileBasics
import SamBasics
import math
import json, zlib, base64
import random

class QualityFormatConverter:
  def __init__(self,type):
    self.observed_to_Q = {}
    self.observed_to_probability = {}
    self.type = type
    if type == 'S':
      start = 33
      for i in range(0,41):
        self.observed_to_Q[i+33] = i
      for i in self.observed_to_Q:
        q = self.observed_to_Q[i]
        self.observed_to_probability[i] = math.pow(10,float(q)/-10)
    elif type == 'P':
      for i in range(0,81):
        self.observed_to_Q[i+33] = i
      for i in self.observed_to_Q:
        q = self.observed_to_Q[i]
        self.observed_to_probability[i] = math.pow(10,float(q)/-10)
    elif type == 'I':
      for i in range(0,41):
        self.observed_to_Q[i+64] = i
      for i in self.observed_to_Q:
        q = self.observed_to_Q[i]
        self.observed_to_probability[i] = math.pow(10,float(q)/-10)
    elif type == 'J':
      for i in range(3,42):
        self.observed_to_Q[i+64] = i
      for i in self.observed_to_Q:
        q = self.observed_to_Q[i]
        self.observed_to_probability[i] = math.pow(10,float(q)/-10)
    elif type == 'L':
      for i in range(0,42):
        self.observed_to_Q[i+33] = i
      for i in self.observed_to_Q:
        q = self.observed_to_Q[i]
        self.observed_to_probability[i] = math.pow(10,float(q)/-10)
    else:
      sys.stderr.write("Error: unsupported quality type "+type+"\n")
      sys.exit()

  def call_observed_ascii_probability(self,ascii_char):
    if ord(ascii_char) not in self.observed_to_probability:
      sys.stderr.write("Error: Looking for a character: '" + ascii_char + "' not present in type: "+self.type+"\n")
      sys.exit()
    return self.observed_to_probability[ord(ascii_char)]

class QualityFormatDetector:
  def __init__(self):
    self.type = 'unknown'
    self.about = 'unknown'
    self.max_read_count = 10000
    self.observed_qualities = {}
  # returns a type that can be used in a quality format converter
  def call_type(self):
    truecount_105_113 = 0
    for i in range(105,114):
      if i in self.observed_qualities: truecount_105_113 += self.observed_qualities[i]
    truecount_76_104 = 0
    for i in range(76,105):
      if i in self.observed_qualities: truecount_76_104 += self.observed_qualities[i]
    truecount_59_63 = 0
    for i in range(59,64):
      if i in self.observed_qualities: truecount_59_63+=self.observed_qualities[i]
    truecount_64_66 = 0
    for i in range(64,67):
      if i in self.observed_qualities: truecount_64_66+=self.observed_qualities[i]
    truecount_67_72 = 0
    for i in range(67,73):
      if i in self.observed_qualities: truecount_67_72+=self.observed_qualities[i]
    truecount_33_58 = 0
    for i in range(33,59):
      if i in self.observed_qualities: truecount_33_58+=self.observed_qualities[i]
    truecount_74 = 0
    if truecount_74 in self.observed_qualities: truecount_74 = self.observed_qualities[74]
    if truecount_74 > 2 and truecount_33_58 >2:
      self.about = "'L' Illumina 1.8+ Phred+33, (0,41) ranges ord 33 to 74"      
      self.type = 'L'
      return 'L'
    if truecount_105_113 > 2 and truecount_33_58 > 2:
      self.about = "'P' PacBio Phred+33, (0,80) ranges ord 33 to 113"      
      self.type = 'P'
      return self.type
    if truecount_33_58 > 2:
      self.about = "'S' Sanger Phred+33, (0,40) ranges ord 33 to 73"      
      self.type = 'S'
      return self.type
    if truecount_59_63 > 2 and truecount_76_104 > 2:
      sys.stderr.write("Warning: Unprogrammed 'X' Solexa Solexa+64, (-5,40) ranges ord 59 to 104\n")      
      return False
    if truecount_64_66 > 2 and truecount_76_104 > 2:
      #sys.stderr.write("Warning: Unprogrammed 'I' Illumina 1.3+ Phred+64, (0,40) ranges ord 64 to 104\n")      
      self.about = "'I' Illumina 1.3+ Phred+64, (0,40) ranges ord 64 to 104"      
      self.type = 'I'
      return self.type
    if truecount_67_72 > 2 and truecount_76_104 > 2:
      #sys.stderr.write("Warning: Unprogrammed 'J' Illumina 1.5+ Phred+64, (3,40) ranges ord 67 to 104\n")      
      self.about = "'J' Illumina 1.5+ Phred+64, (3,40) ranges ord 67 to 104"
      self.type = 'J'
      return self.type
    sys.stderr.write("Warning: unable to choose fastq type\n")
    return False
    

  def set_max_read_count(self,read_count):
    self.max_read_count = read_count

  def read_fastq_file(self,filename):
    gfr = FileBasics.GenericFileReader(filename)
    linecount = 0
    while True and linecount < self.max_read_count:
      line1 = gfr.readline().rstrip()
      if not line1: break
      line2 = gfr.readline().rstrip()
      if not line2: break
      line3 = gfr.readline().rstrip()
      if not line3: break
      line4 = gfr.readline().rstrip()
      if not line4: break
      self.record_observation(line4)
      linecount += 1
    gfr.close()

  # read sam or bam file
  def read_sam_file(self,filename):
    gsr = SamBasics.GenericSamReader(filename)
    linecount = 0
    while True and linecount < self.max_read_count:
      line1 = gsr.readline().rstrip()
      if not line1: break
      line2 = gsr.readline().rstrip()
      if not line2: break
      line3 = gsr.readline().rstrip()
      if not line3: break
      line4 = gsr.readline().rstrip()
      if not line4: break
      self.record_observation(line4)
      linecount += 1
    gsr.close()

  # can be called many times instead of reading a file
  def record_observation(self,line):
    chars = list(line)
    for c in chars: 
      deci = ord(c)
      if deci not in self.observed_qualities:
        self.observed_qualities[deci] = 0
      self.observed_qualities[deci] += 1

# A class to help describe what how qualites are distributed in reads
class QualityProfile:
  def __init__(self,quality_type=None):
    self.quality_type = quality_type
    self.stats = {} #keyed by qual character
    self.observed_lines = 0  # total number of lines looked at
    self.observed_chars = 0  # total number of characters looked at
    self.observed_count_by_position = {} #number times a position was part of an observation
    ### things not needed to be serialized
    self.emitter_tables = None
    self.try_end_runs = True
    return
  def record_observation(self,line):
    chars = list(line.rstrip())
    rlen = len(chars)
    self.observed_lines += 1
    for i in range(0,len(chars)):
      self.observed_chars += 1
      pp = int(100*float(i)/float(rlen)) #position percent
      if pp not in self.observed_count_by_position:
        self.observed_count_by_position[pp] = 0
      self.observed_count_by_position[pp]+=1
      c = chars[i]
      if c not in self.stats:  self.stats[c] = init_stats()
      self.stats[c]['seen'] += 1
      dlen = 0
      diff = 0
      for ci in chars[i:]:
        if ci != c: diff = 1
        if diff == 0: dlen +=1
      # if diff is 0 then it runs to the end
      if diff == 0:
        if pp not in self.stats[c]['runs_to_end_when_seen']:
          self.stats[c]['runs_to_end_when_seen'][pp] = 0
        self.stats[c]['runs_to_end_when_seen'][pp] += 1
      if pp not in self.stats[c]['position']: self.stats[c]['position'][pp] = 0
      self.stats[c]['position'][pp] += 1
      #if diff == 1: #only record length if its not an end run
      if pp not in self.stats[c]['lengths_when_seen']:
        self.stats[c]['lengths_when_seen'][pp] = {}
      if dlen not in self.stats[c]['lengths_when_seen'][pp]:
        self.stats[c]['lengths_when_seen'][pp][dlen] = 0
      self.stats[c]['lengths_when_seen'][pp][dlen] += 1
    #if z % 1000 == 0:
    #  sys.stderr.write(str(z)+"\n")
    #  show_stats(stats)
  def show_stats(self):
    print str(self.quality_type)+ " is quality type"
    print str(self.observed_lines)+" observed QUAL strings"
    print str(self.observed_chars)+" observed QUAL chars"
    #for c in sorted(self.stats.keys()):
    #  print c
    #  for pos in self.stats[c]['position']:
    #    print pos
  def get_serialized(self):
    all = {}
    all['quality_type'] =  self.quality_type
    all['stats'] = self.stats
    all['observed_lines'] = self.observed_lines
    all['observed_chars'] = self.observed_chars
    all['observed_count_by_position'] = self.observed_count_by_position
    str1 = json.dumps(all)
    return encode_name(str1)
  def read_serialized(self,instring):
    all = dict(json.loads(decode_name(instring)))
    self.quality_type = all['quality_type']
    self.stats = all['stats']
    self.observed_lines = all['observed_lines']
    self.observed_chars = all['observed_chars']
    self.observed_count_by_position = {} 
    self.observed_count_by_position = rekey_by_integer(all['observed_count_by_position'])
    for c in self.stats:
      self.stats[c]['runs_to_end_when_seen'] = rekey_by_integer(self.stats[c]['runs_to_end_when_seen'])
      self.stats[c]['position'] = rekey_by_integer(self.stats[c]['position'])
      self.stats[c]['lengths_when_seen'] = rekey_by_integer(self.stats[c]['lengths_when_seen'])
      for pos in self.stats[c]['lengths_when_seen']:
        self.stats[c]['lengths_when_seen'][pos] = rekey_by_integer(	self.stats[c]['lengths_when_seen'][pos])
  def get_probabilities_of_ascii_by_position(self,position):
    #position is a percent (0-99)
    if position not in self.observed_count_by_position:
      position = self.get_nearest_position(position)
    total = self.observed_count_by_position[position]
    #print position
    #print str(total)+" observations at position"
    prob = {}
    for c in sorted(self.stats.keys()):
      if total == 0:
        prob[c] = 0
      if position in self.stats[c]['position']:
        cnt = self.stats[c]['position'][position]
        prob[c] = float(cnt)/float(total)
      else:
        prob[c] = 0
    return prob

  def get_length_distribution_for_ascii_at_position(self,char,position):
    if position not in self.observed_count_by_position:
      position = self.get_nearest_position(position)
    if char not in self.stats: return None
    if position not in self.stats[char]['lengths_when_seen']: return None
    records = 0
    output = {}
    for olen in self.stats[char]['lengths_when_seen'][position]:
      records += self.stats[char]['lengths_when_seen'][position][olen]
    for olen in self.stats[char]['lengths_when_seen'][position]:
      if records == 0: output[olen] = 0
      else:
        output[olen] = float(self.stats[char]['lengths_when_seen'][position][olen])/float(records)
    return output

  def get_end_run_probability_for_ascii_at_position(self,char,position):
    if position not in self.observed_count_by_position:
      position = self.get_nearest_position(position)
    if char not in self.stats: return 0
    if position not in self.stats[char]['position']: return 0
    total = self.stats[char]['position'][position]
    if total == 0: return 0
    if position not in self.stats[char]['runs_to_end_when_seen']: return 0
    endruns = self.stats[char]['runs_to_end_when_seen'][position]
    #print endruns
    #print total
    return float(endruns)/float(total)

  def get_nearest_position(self,position):
    nums = sorted(self.observed_count_by_position.keys())
    nearest = 0
    if position in self.observed_count_by_position: return position
    for i in nums: 
      if i <= position: nearest = i
    return nearest

  #Use all the observed data to make some tables to speed up emitting 
  def initialize_emitter_tables(self):
    self.emitter_tables = {}
    self.emitter_tables['char_by_pos'] = {}
    for pos in self.observed_count_by_position.keys():
      p = self.get_probabilities_of_ascii_by_position(pos)
      tot = 0
      values = []
      for c in sorted(p,key=lambda x:p[x]):
        tot += p[c]
        values.append([tot,c])
      self.emitter_tables['char_by_pos'][pos] = values
    self.emitter_tables['end_run_by_pos_char'] = {}
    for pos in self.observed_count_by_position:
      self.emitter_tables['end_run_by_pos_char'][pos] = {}
      for c in self.stats:
        self.emitter_tables['end_run_by_pos_char'][pos][c] = self.get_end_run_probability_for_ascii_at_position(c,pos)
    self.emitter_tables['run_length_by_pos_char'] = {}
    for pos in self.observed_count_by_position:
      self.emitter_tables['run_length_by_pos_char'][pos] = {}
      for c in self.stats:
        lendist = self.get_length_distribution_for_ascii_at_position(c,pos)
        if lendist:
          tot = 0
          values = []
          for olen in sorted(lendist,key=lambda x:lendist[x]):
            tot += lendist[olen]
            values.append([tot,olen])
          self.emitter_tables['run_length_by_pos_char'][pos][c] = values
      
  def emit(self,rlen):
    if not self.emitter_tables:
      self.initialize_emitter_tables()
    seq = ''
    while len(seq) < rlen:
      seq += self.emit_by_position(len(seq)+1,rlen)
    seq = seq[0:rlen]
    return seq

  #This is an index-1 position, and rlen is the read length
  def emit_by_position(self,pos,rlen):
    if not self.emitter_tables:
      self.initialize_emitter_tables()
    # Get the base
    pp = int(100*float(pos)/float(rlen))
    pp = self.get_nearest_position(pp) # make sure its one we have data for
    rnum = random.random()
    c = None
    #print pp
    #print self.emitter_tables['char_by_pos'].keys()
    for val in self.emitter_tables['char_by_pos'][pp]:
      if rnum < val[0]:
        c = val[1]
        break
    if not c:
      sys.stderr.warning("WARNING: trouble finding quality\n")
    # See if it should be an end run
    #if len(outchars) > 6:
    if self.try_end_runs:
      rnum = random.random()
      if rnum < self.emitter_tables['end_run_by_pos_char'][pp][c]: #an end run
        return c*(rlen-pos+1)
    # See what length it should be
    rnum = random.random()
    if c not in self.emitter_tables['run_length_by_pos_char'][pp]:
      return c
    ldist = self.emitter_tables['run_length_by_pos_char'][pp][c]
    olen = 1
    for val in ldist:
      if rnum < val[0]:
        olen = val[1]
        break
    outchars = c*olen
    #print pp
    #print c
    #print self.emitter_tables['run_length_by_pos_char'][pp][c]
    #print len(outchars)
    #print '---'
    return outchars

def init_stats():
  stats = {}
  stats['seen'] = 0 # counts total
  stats['position'] = {} #counts keyed by position in read
  stats['runs_to_end_when_seen'] = {} #keyed by position in read
  stats['lengths_when_seen'] = {} #keyed by length when seen
  return stats

def encode_name(conversion_string):
  compressed_string = zlib.compress(conversion_string,9)
  enc_string = base64.b32encode(compressed_string)
  return 'QUAL_'+enc_string.rstrip('=')

def decode_name(safename):
  frag = safename.lstrip('QUAL_')
  padding = (8-(len(frag) % 8)) % 8
  c = base64.b32decode(frag+'='*padding)
  return zlib.decompress(c)

def mean(data):
  n = len(data)
  if n < 1:
    sys.stderr.write("Need one or more data point to get a mean\n")
    sys.exit()
  return float(sum(data))/float(n) 

def _ss(data):
  c = mean(data)
  ss = sum([(x-c)**2 for x in data])
  return ss

def stddev(data):
  n = len(data)
  if n < 2:
    sys.stderr.write("Need 2 or more data points for stddev\n")
    sys.exit()
  ss = _ss(data)
  pvar = ss/n 
  return pvar**0.5

def rekey_by_integer(ref):
  temp = {}
  for k1 in ref:
    temp[int(k1)] = ref[k1]
  return temp
