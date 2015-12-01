#!/usr/bin/python
import sys,argparse, re
from subprocess import Popen, PIPE
from SamBasics import MultiEntrySamReader, SAMtoPSLconversionFactory, PSLtoSAMconversionFactory
from PSLBasics import PSL
from SequenceBasics import rc, read_fasta_into_hash

def main():
  parser = argparse.ArgumentParser(description="Take a sam file and join together mate pairs into single alignments.  Alignments must be ordered by query name.")
  parser.add_argument('input',help="FILENAME input .sam or .bam or '-' for STDIN sam")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--sam',action='store_true')
  group.add_argument('--bam',action='store_true')
  parser.add_argument('--mates_only',action='store_true',help="Only output combined mates")
  args = parser.parse_args()
  inf = sys.stdin
  if args.bam or (not args.sam and not args.input == '-'):
    fh = open(args.input)
    p = Popen('samtools view - -h'.split(),stdin=fh,stdout=PIPE)
    inf = p.stdout
  msr = MultiEntrySamReader(inf)
  spc = SAMtoPSLconversionFactory()
  psc = PSLtoSAMconversionFactory()
  # set the headers for the spc
  for h in msr.header:
    print h.rstrip()
    spc.read_header_line(h)
  while True:
    entries = msr.read_entries()
    if not entries: break
    l = []
    r = []
    for sam in entries:
      #Print line if its not a pair
      if not_a_mate_sam(sam):
        if not args.mates_only:
          print sam.get_line()
        continue
      if sam.check_flag(64): l.append(sam)
      if sam.check_flag(128): r.append(sam)
    if not (len(l)==1 and len(r)==1):
      # more than just a unique pair here
      if not args.mates_only:
        for sam in l:  print sam.get_line()
        for sam in r:  print sam.get_line()
      continue
    #Verify pairing by reference and direction
    if l[0].value('rname') != r[0].value('rname') or l[0].check_flag(16) == r[0].check_flag(16):
      sys.stderr.write("ERROR, these are not actually properly paired as we were led to believe\n")
      sys.exit()
    p1 = PSL(spc.convert_line(l[0].get_line()))
    if not re.search('[HP]',l[0].value('cigar')): 
      p1.set_query(l[0].value('seq'))
      p1.set_quality_seq(l[0].value('qual'))
      if l[0].check_flag(16):
        # set the query to what it actually is
        p1.set_query(rc(l[0].value('seq')))
        p1.set_quality_seq(l[0].value('qual')[::-1])      
    p2 = PSL(spc.convert_line(r[0].get_line()))
    if not re.search('[HP]',r[0].value('cigar')): 
      p2.set_query(r[0].value('seq'))
      p2.set_quality_seq(r[0].value('qual'))
      if r[0].check_flag(16):
        # set the query to what it actually is
        p2.set_query(rc(r[0].value('seq')))
        p2.set_quality_seq(r[0].value('qual')[::-1])      
    p12 = join_mated(p1,p2)
    if not p12:
      if not args.mates_only:
        print l[0].get_line()
        print r[0].get_line()
      continue
    #if p1.value('strand') == '-' and p2.value('strand') == '+' \
    #and p2.value('tEnd') < p1.value('tStart'):
    sline = psc.convert_line(p12.get_line(),query_sequence=p12.get_query(),quality_sequence=p12.get_quality_seq())
    #print p12.get_line()
    print sline

def join_mated(p1,p2):
  if p1.value('strand') == '+' and p1.value('tStart') > p2.value('tStart'):
    return False
  if p1.value('strand') == '-' and p1.value('tStart') < p2.value('tStart'):
    return False

  #lets assume p1 is left and p2 is right
  p2r = p2.rc()
  conc = p1.concatonate_queries(p2r)
  return conc

def not_a_mate_sam(sam):
  if not sam.check_flag(1) or not sam.check_flag(2):
    return True
  #see if its not aligned
  if sam.check_flag(4):
    return True
  if sam.check_flag(8):
    return True
        
if __name__=="__main__":
  main()
