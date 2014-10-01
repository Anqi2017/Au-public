import os
import multiprocessing
import re

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
