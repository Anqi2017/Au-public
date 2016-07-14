#!/usr/bin/python
import sys, argparse
from multiprocessing import cpu_count
from subprocess import Popen, PIPE

def main():
  parser = argparse.ArgumentParser(description="Sort multiple bam files, use extension .sorted.bam",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',nargs='+',help="Bam file names")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="thread count")
  args = parser.parse_args()
  for iname in args.input:
    if iname[-4:] != '.bam': 
      sys.stderr.write("ERROR. input must be .bam file\n")
      sys.exit()
  for iname in args.input:
    ofname = iname[:-4]+'.sorted'
    cmd = 'samtools sort -@ '+str(args.threads)+' '+iname+' '+ofname
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split())
    p.communicate()
if __name__=="__main__":
  main()
