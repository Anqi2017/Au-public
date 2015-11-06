#!/usr/bin/python
import argparse, sys, random, os, subprocess
from shutil import rmtree

# Wrapper for a pipeline to convert BWA-mem results into a psl file


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("bwa_bam",help="BAMFILE or - for sam streamed to stdin")
  parser.add_argument("query_sequences",help="FASTA/Q for query sequences")
  parser.add_argument("reference_fasta",help="FASTA for reference genome")
  parser.add_argument("-S","--size",help="FASTA for reference genome")
  group = parser.add_mutually_exclusive_group()
  group.add_argument("--tempdir",default='/tmp',help="DIR to store a temp directory")
  group.add_argument("--specific_tempdir",help="Exact DIR to work in, is not deleted")
  args = parser.parse_args()

  # Manage your tempdir
  # put it into args.tempdir
  if not args.specific_tempdir:
    rnum = random.randint(1,100000000)
    args.tempdir = args.tempdir.rstrip('/')+'/weirathe.'+str(rnum) 
  else:
    args.tempdir = args.specific_tempdir.rstrip('/')
  if not os.path.exists(args.tempdir):
    os.makedirs(args.tempdir)

  # 1. Now we can start the process.  First convert the sam to a psl
  if args.bwa_bam != '-':
    cmd = "samtools view "+args.bwa_bam+" | sam_to_psl.py - | sort_psl.py - --tempdir "+args.tempdir+" -o "+args.tempdir+"/1.psl"
    subprocess.call(cmd,shell=True)
  else:
    cmd = "sam_to_psl.py - | sort_psl.py - --tempdir "+args.tempdir+" -o "+args.tempdir+"/1.psl"
    p1 = subprocess.Popen(cmd,shell=True,stdin=sys.stdin)
    p1.communicate()

  # 2. 

  # Clean up your temporary directory if you aren't in a specific one.
  if not args.specific_tempdir:
    rmtree(args.tempdir)
  return

if __name__=="__main__":
  main()

