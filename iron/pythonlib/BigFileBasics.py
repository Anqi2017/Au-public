import os, sys

# Read lines from a big file starting on the Nth chunk.
# Lets you work on a file as if you had split it, without splitting it.
#
# Pre:  A file name of an Ascii file
#       Needs to have either set_chunk_size_bytes or set_chunk_count run.
#       Example:
#
# Methods:
#       set_chunk_size_bytes(byte_size int)
#       set_chunk_count(chunk_count int)
#       open_fasta_chunk(iter int)
#       open_fastq_chunk(iter int)
#       open_chunk(iter int)
#
# Example:
#       bfcr = BigFileBasics.BigFileChunkReader('my_test_file.txt')
#       bfcr.set_chunk_count(3000)
#       for i in range(0,bfcr.chunk_count):
#         oc = bfcr.open_chunk(i)
#         while True:
#           if not line:
#             break
#           line = oc.read_line()
#           print line.rstrip()
#         oc.close()
#

class BigFileChunkReader:
  def __init__(self,filename):
    self.filename = filename
    self.file_size = os.path.getsize(filename)
    self.chunk_size = self.file_size
    self.chunk_count = 1

  def set_chunk_size_bytes(self,bytes):
    if bytes <= 1:
      sys.stderr.write("Error: chunk_size is too small. "+str(bytes)+" I mean, come on.\n")
      sys.exit()
    self.chunk_size = bytes
    self.chunk_count = int(self.file_size)/int(self.chunk_size)
    mod = int(self.file_size) % int(self.chunk_size)
    if mod != 0: 
      self.chunk_count += 1

  def set_chunk_count(self,chunk_count):
    if chunk_count < 1 or chunk_count >= self.file_size:
      sys.stderr.write("Error: chunk_count is too strange. "+str(chunk_count)+" I mean, come on.\n")
      sys.exit()
    self.chunk_count = chunk_count
    if int(self.file_size) % int(self.chunk_count) == 0:
      self.chunk_size = int(self.file_size)/int(self.chunk_count)
    else:
      self.chunk_size = (int(self.file_size)/int(self.chunk_count))+1

  def open_chunk(self,iter):
    return BigFileChunkReader.OpenChunk(self,iter)

  def open_fasta_chunk(self,iter):
    return BigFileChunkReader.OpenFastaChunk(self,iter)

  def open_fastq_chunk(self,iter):
    return BigFileChunkReader.OpenFastqChunk(self,iter)

  class OpenChunk():
    def __init__(self,parent,iter):
      self.inf = open(parent.filename)
      self.iter = iter
      self.file_size = parent.file_size
      self.chunk_size = parent.chunk_size
      self.has_line_start = False
      self.end_of_chunk = True
      start_pos = self.iter*self.chunk_size
      self.inf.seek(self.iter*self.chunk_size)
      if not at_beginning_of_line(self.inf):
        at_line_start = self.seek_line_start()
      else:
        self.has_line_start = True
      if self.has_line_start == False:
        self.end_of_chunk = True

    def close(self):
      self.inf.close()    

    def read_line(self):
      if self.has_line_start == False:
        return None
      if self.inf.tell() < (self.iter+1)*self.chunk_size:
        return self.inf.readline()
      else:
        self.end_of_chunk = True
        return None

    def seek_line_start(self):
      start_position = self.inf.tell()
      z = 0
      curr = ''
      has_newline = False
      while z < self.chunk_size and self.inf.tell() <= self.file_size and curr != "\n":
        z += 1
        curr = self.inf.read(1)
        if curr == "\n":
          has_newline = True
          self.has_line_start = True
      return

  class OpenFastaChunk():
    def __init__(self,parent,iter):
      self.inf = open(parent.filename)
      self.iter = iter
      self.file_size = parent.file_size
      self.chunk_size = parent.chunk_size
      self.has_line_start = False
      self.end_of_chunk = True
      self.still_in_entry = True
      start_pos = self.iter*self.chunk_size
      self.inf.seek(self.iter*self.chunk_size)
      if not at_beginning_of_fasta_entry(self.inf):
        at_line_start = self.seek_line_start()
      if at_beginning_of_fasta_entry(self.inf):
        self.has_line_start = True
      if self.has_line_start == False:
        self.end_of_chunk = True

    def close(self):
      self.inf.close()    

    def read_line(self):
      if self.has_line_start == False:
        return None
      if self.inf.tell() < (self.iter+1)*self.chunk_size:
        return self.inf.readline()
      # we are past our chunk size
      if self.inf.tell() >= self.file_size:
        self.still_in_entry = False
        self.end_of_chunk = True
        return None
      if at_beginning_of_fasta_entry(self.inf):
        self.still_in_entry = False
        self.end_of_chunk = True
        return None        
      return self.inf.readline()

    def seek_line_start(self):
      start_position = self.inf.tell()
      z = 0
      curr = ''
      has_newline = False
      while z < self.chunk_size and self.inf.tell() <= self.file_size and self.has_line_start == False:
        z += 1
        curr = self.inf.read(1)
        if at_beginning_of_fasta_entry(self.inf):
          has_newline = True
          self.has_line_start = True
      return

  class OpenFastqChunk():
    def __init__(self,parent,iter):
      self.inf = open(parent.filename)
      self.iter = iter
      self.file_size = parent.file_size
      self.chunk_size = parent.chunk_size
      self.has_line_start = False
      self.end_of_chunk = True
      self.still_in_entry = True
      start_pos = self.iter*self.chunk_size
      self.inf.seek(self.iter*self.chunk_size)
      if not at_beginning_of_fastq_entry(self.inf):
        at_line_start = self.seek_line_start()
      if at_beginning_of_fastq_entry(self.inf):
        self.has_line_start = True
      if self.has_line_start == False:
        self.end_of_chunk = True

    def close(self):
      self.inf.close()    

    def read_line(self):
      if self.has_line_start == False:
        return None
      if self.inf.tell() < (self.iter+1)*self.chunk_size:
        return self.inf.readline()
      # we are past our chunk size
      if self.inf.tell() >= self.file_size:
        self.still_in_entry = False
        self.end_of_chunk = True
        return None
      if at_beginning_of_fastq_entry(self.inf):
        self.still_in_entry = False
        self.end_of_chunk = True
        return None        
      return self.inf.readline()

    def seek_line_start(self):
      start_position = self.inf.tell()
      z = 0
      curr = ''
      has_newline = False
      while z < self.chunk_size and self.inf.tell() <= self.file_size and self.has_line_start == False:
        z += 1
        curr = self.inf.read(1)
        if at_beginning_of_fastq_entry(self.inf):
          has_newline = True
          self.has_line_start = True
      return


def at_beginning_of_line(inf):
  if inf.tell() == 0:
    return True
  inf.seek(inf.tell()-1)
  v = inf.read(1)
  if v == "\n":
    return True
  return False

def at_beginning_of_fasta_entry(inf):
  # case for beginning of file
  if inf.tell() == 0:
    next_char = inf.read(1)
    inf.seek(inf.tell()-1)
    if next_char == '>':
      return True
    else:
      return False
  # other cases somewhere in the file
  inf.seek(inf.tell()-1)
  v = inf.read(1)
  if v == "\n":
    next_char = inf.read(1)
    inf.seek(inf.tell()-1)
    if next_char == '>':
      return True
    else:
      return False
  return False

def at_beginning_of_fastq_entry(inf):
  # case for beginning of file
  if inf.tell() == 0:
    next_char = inf.read(1)
    inf.seek(inf.tell()-1)
    if next_char == '@':
      curr_pos = inf.tell()
      l1 = inf.readline()
      l2 = inf.readline()
      second_check = inf.read(1)
      inf.seek(curr_pos)
      if l1 and l2 and second_check == '+':
        return True
      else:
        return False
    else:
      return False
  # other cases somewhere in the file
  inf.seek(inf.tell()-1)
  v = inf.read(1)
  if v == "\n":
    next_char = inf.read(1)
    inf.seek(inf.tell()-1)
    if next_char == '@':
      curr_pos = inf.tell()
      l1 = inf.readline()
      l2 = inf.readline()
      second_check = inf.read(1)
      inf.seek(curr_pos)
      if l1 and l2 and second_check == '+':
        return True
      else:
        return False
  return False

