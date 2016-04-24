#!/usr/bin/python
import sys, argparse, re
from subprocess import Popen, PIPE
from SamBasics import MultiEntrySamReader

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="Bam file in order of query name - for stdin")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-': 
    cmd = "samtools view -h "+args.input
    p = Popen(cmd.split(),stdout=PIPE)
    inf = p.stdout
  mesr = MultiEntrySamReader(inf)
  while True:
    entries = mesr.read_entries()
    if not entries: break
    if len(entries) == 0: break
    if entries[0].value('cigar') == '*': 
      print entries[0].value('qname')+"\t0"
      continue
    sam = entries[0]
    m = re.search('NH:i:(\d+)',sam.entry['remainder'])
    if not m:
      sys.stderr.write("ERROR not a hisat entry\n")
      sys.exit()
    cnt = max([len(entries),int(m.group(1))])
    print entries[0].value('qname')+"\t"+str(cnt)
if __name__=="__main__":
  main()
