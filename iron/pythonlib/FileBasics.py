import re, os, sys
from random import randint
from shutil import rmtree
import subprocess

# A linux stream reader for zipped or unzipped file streams
# Input: A filename
class GenericFileReader:
  def __init__(self,filename):
    self.filename = filename
    self.type  = 'normal'
    if re.search('\.gz$',self.filename): # do it as a gzipped stream
      cmd = 'zcat '+filename
      args = cmd.split()
      self.process = subprocess.Popen(args,stdout=subprocess.PIPE)
      self.type = 'gzipped'
    else:
      self.normal_filehandle = open(filename)
      #cmd = 'cat '+filename
      #args = cmd.split()
      #self.process = subprocess.Popen(args,bufsize=0,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

  def close(self):
    if self.type == 'gzipped':
      if self.process:
        self.process.kill()
    else:
      self.normal_filehandle.close()

  def readline(self):
    if self.type == 'gzipped':
      return self.process.stdout.readline()
    return self.normal_filehandle.readline()

# make_tempdir2
# Makes an empty random directory in tmp and returns the path
# Pre: <basename1> <basename2>
# Post: Creates an empty temporary directory
#       /tmp/<basename1>/<basename2>45923847052/
#       and returns the path as a string
# Modifies: can delete the temporary directory if it already exist

def make_tempdir2(base1,base2):
  if not os.path.isdir("/tmp"):
    print 'file_basics.py error no /tmp'
    sys.exit()
  m = re.match('^[^\/]+$',base1)
  if not m: 
    print 'file_basics.py error bad name base1 '+base1
  m = re.match('^[^\/]+$',base2)
  if not m: 
    print 'file_basics.py error bad name base2 '+base2
  if not os.path.isdir("/tmp/"+base1):
    sys.stderr.write('file_basics.py no temp /tmp/'+base1+" so making it\n")
    os.mkdir("/tmp/"+base1)
  rnum = str(randint(1,100000000))
  tdir = '/tmp/'+base1+'/'+base2+rnum
  if os.path.isdir(tdir):
    sys.stderr.write('file_basics.py strange found '+tdir+"\n")
    sys.stderr.write('file_basics.py should not have... deleting it'+"\n")
    rmtree(tdir)
  os.mkdir(tdir)
  sys.stderr.write('file_basics.py made '+tdir+"\n")
  return tdir
