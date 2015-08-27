#!/usr/bin/python
import argparse, os, sys, re
from subprocess import call

def main():
  parser = argparse.ArgumentParser(description="place a sorted bam file where a sam file is.")
  parser.add_argument('input',help="SAMFILE input file")
  args = parser.parse_args()
  if not os.path.isfile(args.input):
    sys.stderr.write("ERROR: input file not found\n")
    sys.exit()
  m = re.search('(.*)\.sam$',args.input)
  if not m:
    sys.stderr.write("ERROR: wrong input file format\n")
    sys.exit()
  cmd = 'samtools view -Sb '+args.input+' | samtools sort - '+m.group(1)+'.sorted' 
  sys.stderr.write(cmd+"\n")
  call(cmd,shell=True)

if __name__=="__main__":
  main()
