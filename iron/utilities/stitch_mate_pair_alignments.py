#!/usr/bin/python
import sys,argparse, re
from subprocess import Popen, PIPE
from SamBasics import MultiEntrySamReader, SAMtoPSLconversionFactory, PSLtoSAMconversionFactory
from PSLBasics import PSL
from SequenceBasics import rc, read_fasta_into_hash
from multiprocessing import Pool, cpu_count, Lock

#Global
glock = Lock()

def main():
  parser = argparse.ArgumentParser(description="Take a sam file and join together mate pairs into single alignments.  Alignments must be ordered by query name.")
  parser.add_argument('input',help="FILENAME input .sam or .bam or '-' for STDIN sam")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--sam',action='store_true')
  group.add_argument('--bam',action='store_true')
  parser.add_argument('--mates_only',action='store_true',help="Only output combined mates")
  parser.add_argument('--threads',type=int,default=1,help="Number of threads to use, default is 1")
  args = parser.parse_args()
  inf = sys.stdin
  if args.bam or (not args.sam and not args.input == '-'):
    fh = open(args.input)
    p = Popen('samtools view - -h'.split(),stdin=fh,stdout=PIPE)
    inf = p.stdout
  buffer_size = 10000
  buffer = []
  msr = MultiEntrySamReader(inf)
  spc = SAMtoPSLconversionFactory()
  psc = PSLtoSAMconversionFactory()
  # set the headers for the spc
  for h in msr.header:
    print h.rstrip()
    spc.read_header_line(h)
  if args.threads > 1:
    p1 = Pool(processes=args.threads)
  while True:
    entries = msr.read_entries()
    if not entries: break
    buffer.append(entries)
    if len(buffer) >= buffer_size:
      if args.threads > 1:
        p1.apply_async(do_buffer,args=(buffer,msr,spc,psc,args),callback=do_callback)
      else: 
        v = do_buffer(buffer,msr,spc,psc,args)
        do_callback(v)
      buffer = []
  if len(buffer) > 0:
    if args.threads > 1:
      p1.apply_async(do_buffer,args=(buffer,msr,spc,psc,args),callback=do_callback)
    else:
      v = do_buffer(buffer,msr,spc,psc,args)
      do_callback(v)
  if args.threads > 1:
    p1.close()
    p1.join()

def do_callback(outputs):
  global glock
  glock.acquire()
  for output in outputs:
    print output
  glock.release()

def do_buffer(buffer,msr,spc,psc,args):
  outputs = []
  for entries in buffer:
    l = []
    r = []
    for sam in entries:
      #Print line if its not a pair
      if not_a_mate_sam(sam):
        if not args.mates_only:
          outputs.append(sam.get_line())
        continue
      if sam.check_flag(64): l.append(sam)
      if sam.check_flag(128): r.append(sam)
    if not (len(l)==1 and len(r)==1):
      # more than just a unique pair here
      if not args.mates_only:
        for sam in l:  outputs.append(sam.get_line())
        for sam in r:  outputs.append(sam.get_line())
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
        outputs.append(l[0].get_line())
        outputs.append(r[0].get_line())
      continue
    #if p1.value('strand') == '-' and p2.value('strand') == '+' \
    #and p2.value('tEnd') < p1.value('tStart'):
    sline = psc.convert_line(p12.get_line(),query_sequence=p12.get_query(),quality_sequence=p12.get_quality_seq())
    #print p12.get_line()
    outputs.append(sline)
  return outputs

def join_mated(p1,p2):
  if p1.value('strand') == '+' and p1.value('tStart') > p2.value('tStart'):
    return False
  if p1.value('strand') == '-' and p1.value('tStart') < p2.value('tStart'):
    return False

  #lets order them here 
  left = p1
  right = p2
  if p1.value('tStart') > p2.value('tStart'):
    left = p2
    right = p1
  rightrc = right.rc()
  conc = left.concatonate_queries(rightrc)
  #print conc
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
