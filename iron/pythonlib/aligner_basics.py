import os, sys
import multiprocessing
import re

import file_basics
from shutil import rmtree, copyfile

# pre:       A fasta file name, an output base name
# post:      Index is built at the base name
# modifies:  Creates index files

def build_gmap_index(fasta_filename,base_name):
  dir = os.path.dirname(base_name)
  m = re.search('([^\/]+)$',base_name)
  if not m:
    print 'no valid base name for build_gmap_index '+basename
    sys.exit()
  name = m.group(1)
  cmd = 'gmap_build -D '+dir+' -d '+name+' ' + fasta_filename
  os.system(cmd)

# pre: <reads filename> <bowtie2 index base name> <samfile output name>
#      If the read name is a .fa or .fasta it will switch to fasta mode
#      otherwise it expects a fastq format file.
# post: writes a sam file out
# modifies: reads multiprocessing to choose cpu count. 

def gmap_all(reads_filename,index_base_name,outpsl_filename):
  tdir = file_basics.make_tempdir2('weirathe','gmap')
  readformattag = ''
  corecount = str(multiprocessing.cpu_count())
  m = re.match('^(.*)\/([^\/]+)$',index_base_name)
  if not m:
    print "error: path should include both a directory and basename you should be able to use ./mybase if its in the same directory you are currently in"
    sys.exit()
  cmd = 'gmap -D '+m.group(1)+' -f 1 -d '+m.group(2)+' -t '+corecount+' '+reads_filename+' 1> '+tdir+'/all.psl 2>/dev/null'
  sys.stderr.write(cmd+"\n")
  os.system(cmd)  
  copyfile(tdir+'/all.psl',outpsl_filename)
  rmtree(tdir)

# pre:       A fasta file name, an output base name
# post:      Index is built at the base name
# modifies:  Creates index files

def build_bowtie2_index(fasta_filename,base_name):
  cmd = 'bowtie2-build '+fasta_filename+' '+base_name
  os.system(cmd)

# pre: <reads filename> <bowtie2 index base name> <samfile output name>
#      If the read name is a .fa or .fasta it will switch to fasta mode
#      otherwise it expects a fastq format file.
# post: writes a sam file out
# modifies: reads multiprocessing to choose cpu count. 

def bowtie2_unpaired(reads_filename,index_base_name,sam_filename):
  readformattag = ''
  if re.search('\.fa$|\.fasta',reads_filename):
    readformattag = ' -f '
  corecount = multiprocessing.cpu_count()
  cmd = 'bowtie2 -k 10 -p '+str(corecount)+readformattag+' -x '+index_base_name+' -U '+reads_filename+' > '+sam_filename
  os.system(cmd)  

def get_make_tempdir():
  if not os.path.isdir("/tmp/"):
    print 'no tmp'
    sys.exit()
  if not os.path.isdir("/tmp/weirathe"):
    print 'no temp weirathe'
    os.mkdir("/tmp/weirathe")
  rnum = str(randint(1,100000000))
  #remove this next line later "rnum = str(1000)"
  rnum = str(1000)
  tdir = '/tmp/weirathe/anlong'+rnum
  #put back later
  #if os.path.isdir(tdir):
  #  print 'strange found '+tdir
  #  print 'should not have... deleting it'
  #  rmtree(tdir)
  #os.mkdir(tdir)
  return tdir
