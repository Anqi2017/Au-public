#!/usr/bin/python
import sys, argparse, re, gzip
from subprocess import Popen, PIPE
from Bio.Format.Sam import SAM

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use BAM or - for SAM STDIN")
  parser.add_argument('-o','--output',help="output file")
  args = parser.parse_args()
  
  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')

  inf = sys.stdin
  if args.input != '-':
    cmd = "samtools view "+args.input
    p = Popen(cmd.split(),stdout=PIPE)
    inf = p.stdout
  for line in inf:
    v = SAM(line)
    for rng in [x.get_range() for x in v.get_target_transcript(68).exons]:
      of.write("\t".join([str(x) for x in rng.get_bed_array()])+"\n")

  if args.input != '-':
    p.communicate()    

  of.close()

if __name__=="__main__":
  main()
