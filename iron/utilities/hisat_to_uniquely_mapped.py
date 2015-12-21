#!/usr/bin/python
import sys, argparse, re
from subprocess import Popen, PIPE
from SamBasics import SAM, is_header

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="BAM FILE or '-' for STDIN SAM format")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-': 
    cmd = "samtools view -F 4 -h "+args.input
    p = Popen(cmd.split(),stdout=PIPE)
    inf = p.stdout
  for line in inf:
    if is_header(line): 
      print line.rstrip()
      continue
    sam = SAM(line)
    if sam.entry['cigar'] == '*': continue
    m = re.search('NH:i:(\d+)',sam.entry['remainder'])
    if not m:
      sys.stderr.write("ERROR not a hisat entry\n")
      sys.exit()
    if int(m.group(1))==1:
      print line.rstrip()
if __name__=="__main__":
  main()
