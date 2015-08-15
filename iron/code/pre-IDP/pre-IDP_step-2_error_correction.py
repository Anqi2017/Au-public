#!/usr/bin/python
import argparse, sys, os, random, multiprocessing, subprocess, re
from shutil import rmtree, copytree

# This will operate a pipeline to go from 
# Pre:
#     pacbio_raw .. the path to the bax.h5 or bas.h5 you want to process
#                please keep in mind that in this requires the accordingly 
#                named xml file to be in the parent directory as the results
#                are stored by pacbio.
#     output directory must not already exist.
# Post:
#     output direcotry by default is in the current working directory
#     and called pre-IDP_output/ but can be anything by settting --output
#         

def main():
  parser = argparse.ArgumentParser(description="take the long and short reads to perform LSC",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--short_reads',nargs='+',required=True,help='READFILE or READFILE_1 READFILE_2 for paired end')
  parser.add_argument('--short_read_type',required=True,choices=['FASTA','FASTQ'],help='FASTA or FASTQ')
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--step_1_folder',help='FOLDERNAME of step one outputs containing long reads to correct')
  group.add_argument('--long_read_fasta_to_correct',help='FASTAFILE input long reads')
  parser.add_argument('--threads',type=int,default=0,help='INT number of threads to use')
  group1 = parser.add_mutually_exclusive_group()
  group1.add_argument('--tempdir',default='/tmp',help='FOLDERNAME location of random temporary directory')
  group1.add_argument('--specific_tempdir',help='FOLDERNAME location of the specific temporary directory.  will not be removed on completion.')
  parser.add_argument('--output',default='output_pre-IDP_step-2_error_correction',help='FOLDERNAME of output must not already exist')
  parser.add_argument('--save_tempdir',help='DIRECTORYNAME name of a directory to be created that has the temporary folder in it.')
  args = parser.parse_args()

  if args.threads==0:
    args.threads = multiprocessing.cpu_count()

  if args.short_reads and not args.short_read_type:
    sys.stderr.write("ERROR please specify --short_read_type as FASTA or FASTQ")
    sys.exit()

  if os.path.isdir(args.output):
    sys.stderr.write("ERROR "+args.output+" output folder must not already exist.\n")
    sys.exit()

  if args.save_tempdir:
    if os.path.isdir(args.save_tempdir):
      sys.stderr.write("ERROR "+args.save_tempdir+" folder must not already exist.\n")
      sys.exit()

  tdir = setup_temporary_directory(args)  
  sys.stderr.write("working in "+tdir+"\n")

  long_reads = args.long_read_fasta_to_correct
  if args.step_1_folder:
    long_reads = args.step_1_folder.rstrip('/')+'/ccs_lq_and_longest_subreads_to_correct.fa'

  # now we can perform LSC if we have short reads
  sys.stderr.write("Perform LSC correction on reads\n")
  sr_fasta = tdir+'/SR.fasta'
  of = open(sr_fasta,'w')
  for fname in args.short_reads:
    sys.stderr.write("Fetching "+fname+"\n")
    if args.short_read_type == 'FASTA':
      with open(fname) as inf:
        for line in inf:
          of.write(line)
    elif args.short_read_type == 'FASTQ':
      gfr = GenericFastqFileReader(fname)
      while True:
        entry = gfr.read_entry()
        if not entry: break
        of.write(">"+entry['name']+"\n")
        of.write(entry['seq']+"\n")
    else:
      sys.stderr.write("ERROR unrecongnized short read type\n")
      sys.exit()
  of.close()
  sys.stderr.write("Executing LSC\n")
  
  execute_lsc(tdir,long_reads,sr_fasta,args)
  sys.stderr.write("Finished LSC\n")
  copytree(tdir+'/output',args.output)
  if args.save_tempdir:
    copytree(tdir,args.save_tempdir)
  if not args.specific_tempdir:
    rmtree(tdir)  

def execute_lsc(tdir,to_correct_fasta,sr_fasta,args):
  config_string = '''##
###################################################
#
# This cofiguration file contains all settings for a run
# of LScorr.
#
# lines begining with '#' are comments
# lists begin with '> tag' and end with '<' on separate lines
#
###################################################
##

#########################
## Required Settings
##

##
# python path 
# (single value)

python_path = /usr/bin/python

##
# Run mode 
# (single value)
# 0: end-to-end
# 1: alignment stage (generating compressed SRs to LRs alignments)
# 2: correcting stage (assuming stage 1 is already done )

mode = 0

##
# Long reads file
# (single value)

LR_pathfilename = '''+to_correct_fasta+'''

##
# Long reads file type
# (single value:  fa or fq)

LR_filetype = fa


##
# Short reads file
# (single value)

SR_pathfilename = '''+sr_fasta+'''

##
# Short reads file type
# (single value:  fa or fq or cps)
# If you have run LSC on the same SR data before, you can find the compressed SR data in temp folder (SR.fa.cps and SR.fa.idx files). 
# You can point the SR_pathfilename to SR.fa.cps file (the same folderpath should also include SR.fa.idx file) 
# In this case generating compressed short reads would be skipped
# (single value)

SR_filetype = fa


##
# Is this nonredundant SR data set? (Y or N)
# If you have run LSC on the same SR data before, you could find it in temp folder. Its name is "SR_uniq.fa".
# You can use this "SR_uniq.fa" as the short reads file and set this option to "Y"
# (single value)

I_nonredundant = N

## 
# Short-reads coverage depth (SCD)
# Generates LR-SR alignemnt file with expected SR coverage depth of SCD value.
# Note: SCD filter would be applied to LR segments with SR coverage of more than SCD value. 
# -1: All alignemnt results are used in correction step (no filtration).
# positive integer: Filters SR-LR alignment results to have expected SR coverage depth of SCD. 
# (positive integer or -1)

SCD = 20

##
# Number of threading for short reads alignment to long reads
# (single value)

Nthread1 = '''+str(args.threads)+'''

##
# Number of threading for corrections
# (single value)

Nthread2 = '''+str(args.threads)+'''

##
# Max memory usage for unix sort command (-S option) per thread depending on your system memory limit
# Note: This setting is per thread and number of threads is set through Nthread1 and Nthread2 parameters
# -1: no-setting (default sort option) 
# example: 4G , 100M , ...

sort_max_mem = -1

#########################

##
# Temp folder
# (single value)

temp_foldername = '''+tdir+'/lsc_temp'+'''

##
# Output folder
# (single value)

output_foldername = '''+tdir+'/output'+'''


##
# Remove PacBio tails sub reads? (Y or N)
# The names of PacBio long reads must be in the format of the following example: ">m111006_202713_42141_c100202382555500000315044810141104_s1_p0/16/3441_3479"
# The last two numbers (3441 and 3479 in this example) are the positions of the sub reads. 
# (single value)

I_RemoveBothTails = Y

##
# Min. number of non'N' after compressing 
# (single value)

MinNumberofNonN = 10

##
# Max 'N' are allowed after compressing
# (single value)

MaxN = 2

##
# Remove intermediate  files at the end of LSC run (for instance:  aligner sam output, LR-SR mapping files, ...)
# (0: No, 1: Yes )

clean_up = 1


#########################

##
# Aligner could be set to novoalign, bwa or bowtie2

aligner = bowtie2

# Maximum error rate percentage to filter compressed LR-SR alignment results
# (single value)

max_error_rate = 10

# Aligner command options   
# Note: Do not specify number of threads in the command options, it is set through Nthread1

bowtie2_options = --end-to-end -a -f -L 15 --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.08 --no-unal
bwa_options =  -n 0.08 -o 20 -e 3 -d 0 -i 0 -M 1 -O 0 -E  1 -N 
novoalign_options =  -t 0,1.5 -g 0 -x 20 -r Ex 1000 -R 500 -o Sa 
razers3_options = -i 92'''
  of = open(tdir+'/run.cfg','w')
  of.write(config_string+"\n")
  of.close()
  mypath = os.path.realpath(__file__)
  mydir = '.'
  m = re.match('^(.*)\/[^\/]+\.py',mypath)
  if m: mydir = m.group(1)
  cmd1 = mydir+'/LSC_1_beta/bin/runLSC.py '+tdir+'/run.cfg'
  subprocess.call(cmd1,shell=True)

def setup_temporary_directory(args):
  if args.specific_tempdir:
    tdir = args.specific_tempdir.rstrip('/')
    if not os.path.isdir(tdir):
      os.makedirs(tdir)
    return tdir
  if not os.path.isdir(args.tempdir):
    sys.stderr.write("ERROR invalid temporary directory "+args.tempdir+"\n")
    sys.exit()
  tdir = args.tempdir.rstrip('/')+'/'+'preidp.'+str(random.randint(1,10000000))
  os.makedirs(tdir)
  if not os.path.isdir(tdir):
    sys.stderr.write("ERROR failed to make working temp directory "+args.tempdir+"\n")
    sys.exit()
  return tdir

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


if __name__=="__main__":
  main()
