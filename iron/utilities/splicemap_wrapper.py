#!/usr/bin/python
import os, subprocess, argparse, sys, re
from random import randint
from shutil import rmtree, copytree

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--bowtie_index',nargs=1,help='path to bowtie index',required=True)
  parser.add_argument('--genome',nargs=1,help='path reference genome',required=True)
  parser.add_argument('--read_mismatches',nargs=1,help='Number of read mismatches (default 2)')
  parser.add_argument('--threads',nargs=1,help='Number of threads (default 1)')
  parser.add_argument('--tempdir',nargs=1,help='DIRECTORY location to store temp files')
  parser.add_argument('--output',help='DIRECTORYNAME path of directory to save result')
  parser.add_argument('--output_all',help='DIRECTORYNAME path of directory to save result')
  parser.add_argument('--read_type',default='FASTQ',help='Read type FASTQ or FASTA')
  parser.add_argument('reads',nargs='+',help='reads (second file is for a mate')
  args = parser.parse_args()
  threads = 2
  if args.threads:
    threads = int(args.threads[0])
  read_mismatches = 2
  if args.read_mismatches:
    read_mismatches = int(args.read_mismatches[0])
  if len(args.reads) > 2:
    sys.stderr.write("Too many read files.  Takes only one fastq file or two (a mate pair)\n")
    return
  reads1 = args.reads[0]
  reads2 = None
  if len(args.reads) == 2:
    reads2 = args.reads[1]
  bdir = args.bowtie_index[0]
  #'/Shared/Au/jason/Reference/UCSC/Human/hg19_GRCh37_feb2009/Bowtie_Index/genome'
  gdir = './'
  wcard = args.genome[0]
  m = re.match('^(.*\/)([^\/]+)$',args.genome[0])
  if m:
    gdir = m.group(1)
    wcard = m.group(2)
  #'/Shared/Au/jason/Reference/UCSC/Human/hg19_GRCh37_feb2009/Genome/'
  tstart = '/tmp'
  if args.tempdir:
    tstart = args.tempdir[0]
  tdir = tstart.rstrip('/')+'/'+'weirathe.'+str(randint(1,100000000))
  if not os.path.exists(tdir): os.makedirs(tdir)
  # Make a new reads 1 if its gzipped
  if re.search('.gz$',reads1):
    subprocess.call('zcat '+reads1+' > '+tdir+'/reads1.fq',shell=True)
    reads1 = tdir+'/reads1.fq'
  if re.search('.gz$',reads2):
    subprocess.call('zcat '+reads2+' > '+tdir+'/reads2.fq',shell=True)
    reads2 = tdir+'/reads2.fq'
  cfg = get_cfg(reads1,reads2,bdir,gdir,wcard,threads,read_mismatches,tdir,args.read_type)
  of = open(tdir+'/run.cfg','w')
  of.write(cfg)
  of.close()
  cmd = "runSpliceMap "+tdir+'/run.cfg'
  FNULL = open(os.devnull,'w')
  stream = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE,stderr=FNULL)
  while True:
    e = stream.stdout.readline()
    if not e: break
  if not args.output:
    with open(tdir+'/output/good_hits.sam') as inf:
      for line in inf:
        print line.rstrip()
  if args.output_all:
    copytree(tdir,args.output)
  elif args.output:
    copytree(tdir+'/output/',args.output)
  rmtree(tdir)

def get_cfg(reads1,reads2,bdir,gdir,wcard,threads,read_mismatches,tdir,read_type):
  cfg = '''\
genome_dir = '''+gdir+'''
> reads_list1
'''+reads1+'''
<
'''
  if reads2:
    cfg += '''\
> reads_list2
'''+reads2+'''
<
'''
  cfg += '''\
read_format = '''+read_type+'''
mapper = bowtie
temp_path = '''+tdir+'/temp'+'''
out_path = '''+tdir+'/output'+'''
max_intron = 400000
min_intron = 20000
max_multi_hit = 10
seed_mismatch = 1
read_mismatch = '''+str(read_mismatches)+'''
sam_file = cuff
ud_coverage = yes
chromosome_wildcard = '''+wcard+'''
num_chromosome_together = '''+str(threads)+'''
bowtie_base_dir = '''+bdir+'''
num_threads = '''+str(threads)+'''
try_hard = yes'''
  return cfg


main()
