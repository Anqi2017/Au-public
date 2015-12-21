#!/usr/bin/python
import argparse, sys, os, random, subprocess
from shutil import rmtree

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("fastq_file",help="FASTQFILE or - for stdin")
  parser.add_argument("-o","--output",help="OUTFILE or STDOUT if not specified")
  parser.add_argument("--tempdir",help="DIR to store a temp directory")
  parser.add_argument('-S',"--size",help="S option for linux sort units are kb unless specified")
  args = parser.parse_args()

  # Manage your tempdir
  # put it into args.tempdir
  if args.tempdir:
    rnum = random.randint(1,100000000)
    args.tempdir = args.tempdir.rstrip('/')+'/weirathe.'+str(rnum) 
    if not os.path.exists(args.tempdir):
      os.makedirs(args.tempdir)
  
  inf = sys.stdin
  if args.fastq_file != '-':
    inf = open(args.fastq_file)
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')

  cmd = "fastq_to_tsv.pl | sort -k 1,1"
  if args.tempdir:
    cmd += " -T "+args.tempdir
  if args.size:
    cmd += " -S "+args.size
  cmd += " | tsv_to_fastq.pl"
  p1 = subprocess.Popen(cmd,shell=True,stdin = subprocess.PIPE, stdout = of)
  for line in inf:
    p1.stdin.write(line)
  p1.communicate()
  inf.close()
  of.close()
  #cleanup
  if args.tempdir:
    rmtree(args.tempdir)

if __name__=="__main__":
  main()
