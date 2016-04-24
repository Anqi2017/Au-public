#!/usr/bin/python
import sys, argparse, re, os
from SamBasics import SamLocusStream
from subprocess import Popen, PIPE
from multiprocessing import Pool, Lock

glock = Lock()

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Sorted bam (preferrably indexed) Use - for STDIN sam. If streaming in be sure to remove unmapped reads")
  parser.add_argument('--threads',type=int,default=1,help="use multiple threads the bam has been indexed.  Order is not preserved.")
  args = parser.parse_args()


  single_thread = True
  if args.threads == 1: single_thread = True
  elif args.input != '-':
    if os.path.isfile(args.input+'.bai'): single_thread = False
    else: 
      single_thread = True
      sys.stderr.write("Warning doing single thread because lacking index\n")

  chrs = None
  if args.input != '-':
    chrs = set()
    cmd = 'samtools view -H '+args.input
    p = Popen(cmd.split(),stdout=PIPE)
    for line in p.stdout:
      m = re.match('@SQ\s+SN:(\S+)\s+LN:\d+',line)
      if m: chrs.add(m.group(1))
    p.communicate()
  #easy case of single thread
  if single_thread:
    if args.input == '-':
      dostream(sys.stdin)
    else:
      cmd = 'samtools view -F 4 -h '+args.input
      p = Popen(cmd.split(),stdout=PIPE)
      dostream(p.stdout)
      p.communicate()
  else:
    p = Pool(processes=args.threads)
    for chr in sorted(chrs):
      p.apply_async(dofilestream,args=(args.input,chr),callback=printres)
    p.close()
    p.join()

def printres(res):
  global glock
  glock.acquire()
  for line in res:
    print line
  glock.release()

def dofilestream(filename,chrom):
  cmdx = 'samtools view -h -F 4 '+filename+' '+chrom
  px = Popen(cmdx.split(),stdout=PIPE)
  res = dopart(px.stdout)
  px.communicate()
  return res

def dopart(stream):
    reader = SamLocusStream(stream)
    v = []
    while True:
      locus = reader.read_locus()
      if not locus: break
      v.append(locus[0].chr+"\t"+str(locus[0].start-1)+"\t"+str(locus[0].end))
    return v

def dostream(stream):
    reader = SamLocusStream(stream)
    while True:
      locus = reader.read_locus()
      if not locus: break
      print locus[0].chr+"\t"+str(locus[0].start-1)+"\t"+str(locus[0].end)


if __name__=="__main__":
  main()
