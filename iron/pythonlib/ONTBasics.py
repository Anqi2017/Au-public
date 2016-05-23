# ONBasics.py
#
# Description:
#  Hooks oxford nanopore data
#
import h5py, re, os, sys

class fast5:
  # init
  # Pre: db database ie hg19
  def __init__(self,filename):
    if not os.path.isfile(filename):
      sys.stderr.write("ERROR bad filename for h5: "+filename+"\n")
      sys.exit()
    self.trim_header = False
    self.h5 = h5py.File(filename,'r')
    #self.read_2d = self.extract_2D()
    #self.read_template = self.extract_template()
    #self.read_complement = self.extract_complement()

  def set_trim_header_bool(self,choice):
    self.trim_header = choice

  def close(self):
    self.h5.close()
    
  def extract_2D(self):
    if not self.has_basic_read_info(): return None
    if not 'Basecall_2D_000' in self.h5['Analyses']: return None
    if not 'BaseCalled_2D' in self.h5['Analyses']['Basecall_2D_000']: return None
    return fast5.fastq(self.h5['Analyses']['Basecall_2D_000']['BaseCalled_2D']['Fastq'][()],self.trim_header)
  def extract_complement(self):
    if self.has_basic_R7_read_info():
      if not 'BaseCalled_complement' in self.h5['Analyses']['Basecall_2D_000']: return None
      return fast5.fastq(self.h5['Analyses']['Basecall_2D_000']['BaseCalled_complement']['Fastq'][()],self.trim_header)
    elif self.has_basic_R9_read_info():
      if not 'BaseCalled_complement' in self.h5['Analyses']['Basecall_1D_000']: return None
      return fast5.fastq(self.h5['Analyses']['Basecall_1D_000']['BaseCalled_complement']['Fastq'][()],self.trim_header)
    else:
      return None
  def extract_template(self):
    if self.has_basic_R7_read_info():
      if not 'BaseCalled_template' in self.h5['Analyses']['Basecall_2D_000']: return None
      return fast5.fastq(self.h5['Analyses']['Basecall_2D_000']['BaseCalled_template']['Fastq'][()],self.trim_header)
    elif self.has_basic_R9_read_info():
      if not 'BaseCalled_template' in self.h5['Analyses']['Basecall_1D_000']: return None
      return fast5.fastq(self.h5['Analyses']['Basecall_1D_000']['BaseCalled_template']['Fastq'][()],self.trim_header)
    else:
      return None

  def has_basic_read_info(self):
    if self.has_basic_R9_read_info(): return True
    if self.has_basic_R7_read_info(): return True
    return False

  def has_basic_R9_read_info(self):
    if not 'Analyses' in self.h5: return False
    #if not 'Basecall_2D_000' in self.h5['Analyses']: return False
    if not 'Basecall_1D_000' in self.h5['Analyses']: return False
    #if not 'BaseCalled_complement' in self.h5['Analyses']['Basecall_1D_000']: return False
    if not 'BaseCalled_template' in self.h5['Analyses']['Basecall_1D_000']: return False
    return True
  def has_basic_R7_read_info(self):
    if not 'Analyses' in self.h5: return False
    if not 'Basecall_2D_000' in self.h5['Analyses']: return False
    if not 'BaseCalled_complement' in self.h5['Analyses']['Basecall_2D_000']: return False
    if not 'BaseCalled_template' in self.h5['Analyses']['Basecall_2D_000']: return False
    return True

  class fastq:
    def __init__(self,input,trim_header):
      [name,seq,third,qual] = input.rstrip().split("\n")
      m = re.match('^@(.*)$',name)
      if trim_header:
        m = re.match('^@(\S+)',name)
      self.name = m.group(1)
      self.seq = seq
      self.qual = qual
      return      
    def length(self):
      return len(self.seq)
    def fastq(self):
      return "@"+self.name+"\n"+self.seq +"\n+\n" + self.qual + "\n"
