#!/usr/bin/python
import sys,argparse, re
from SamBasics import MultiEntrySamReader, SAMtoPSLconversionFactory
from PSLBasics import PSL

def main():
  parser = argparse.ArgumentParser(description="Take a sam or bam file and output the best alignment for each read, it still can output the same read name twice if they happen to be mate pairs, but it will only output the best alignment for each individual mate, not necessarily the two together.  You could combine mates if that is helpful with another script.")
  parser.add_argument('input',help="FILENAME input .sam or .bam or '-' for STDIN sam")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--bam',action='store_true')
  group.add_argument('--sam',action='store_true')
  args = parser.parse_args()
  inf = sys.stdin
  if args.bam or (not args.sam and not args.input == '-'):
    fh = open(args.input)
    p = Popen('samtools view - -h'.split(),stdin=fh,stdout=PIPE)
    inf = p.stdout
  msr = MultiEntrySamReader(inf)
  spc = SAMtoPSLconversionFactory()
  # set the headers for the spc
  for h in msr.header:
    print h.rstrip()
    spc.read_header_line(h)
  while True:
    entries = msr.read_entries()
    if not entries: break
    longest0 = 0
    entry0 = None
    longest1 = 0
    entry1 = None
    longest2 = 0
    entry2 = None
    for sam in entries:
      pline = spc.convert_line(sam.get_line())
      if not pline: continue
      side = None
      if sam.check_flag(64): side = 1
      if sam.check_flag(128): side = 2
      p = PSL(pline)
      if p.get_coverage() > longest0:
        longest0 = p.get_coverage()
        entry0 = sam
      if side == 1 and p.get_coverage() > longest1:
        longest1 = p.get_coverage()
        entry1 = sam
      if side == 2 and p.get_coverage() > longest2:
        longest2 = p.get_coverage()
        entry2 = sam
    if entry0: #output the combined if its there
      print entry0.get_line()
    else:
      if entry1: #output each of the mates if they are paired but not joined
        print entry1.get_line()
      if entry2:
        print entry2.get_line()
if __name__=="__main__":
  main()
