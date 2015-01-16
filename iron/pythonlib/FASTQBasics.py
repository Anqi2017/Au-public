import re, sys
import subprocess
import FileBasics
import math

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
      sys.stderr.write("Error: Unprogrammed 'L' Illumina 1.8+ Phred+33, (0,41) ranges ord 33 to 74\n")      
      sys.exit()
      return
    if truecount_33_58 > 2:
      self.about = "'S' Sanger Phred+33, (0,40) ranges ord 33 to 73"      
      self.type = 'S'
      return self.type
    if truecount_59_63 > 2 and truecount_76_104 > 2:
      sys.stderr.write("Error: Unprogrammed 'X' Solexa Solexa+64, (-5,40) ranges ord 59 to 104\n")      
      sys.exit()
      return
    if truecount_64_66 > 2 and truecount_76_104 > 2:
      sys.stderr.write("Error: Unprogrammed 'I' Illumina 1.3+ Phred+64, (0,40) ranges ord 64 to 104\n")      
      sys.exit()
      return
    if truecount_67_72 > 2 and truecount_76_104 > 2:
      print 'J'
      sys.stderr.write("Error: Unprogrammed 'J' Illumina 1.5+ Phred+64, (3,40) ranges ord 67 to 104\n")      
      sys.exit()
      return
    sys.stderr.write("Error: unable to choose fastq type\n")
    sys.exit()
    

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

  def record_observation(self,line):
    chars = list(line)
    for c in chars: 
      deci = ord(c)
      if deci not in self.observed_qualities:
        self.observed_qualities[deci] = 0
      self.observed_qualities[deci] += 1
