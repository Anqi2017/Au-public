#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree
from multiprocessing import cpu_count, Pool, Lock
from subprocess import PIPE, Popen
from tempfile import mkdtemp, gettempdir

glock = Lock()
of = sys.stdout

def main():
  #do our inputs
  args = do_inputs()
  global of
  if args.output:
    of = open(args.output,'w')
  if args.threads > 1:
    p = Pool(args.threads)
  for fname in args.inputs:
    if args.threads > 1:
      p.apply_async(do_bam,args=(fname,),callback=do_out)
    else:
      res = do_bam(fname)
      do_out(res)
  of.close()
  if args.threads > 1:
    p.close()
    p.join()

def do_out(res):
  global glock
  global of
  glock.acquire()
  of.write(res[0]+"\t"+str(res[1])+"\n")
  glock.release()

def do_bam(fname):
  cmd1 = "samtools view "+fname
  cmd2 = "cut -f 1"
  cmd3 = "sort"
  cmd4 = "uniq"
  cmd5 = "wc -l"
  p1 = Popen(cmd1.split(),stdout=PIPE)
  p2 = Popen(cmd2.split(),stdin=p1.stdout,stdout=PIPE)
  p3 = Popen(cmd3.split(),stdin=p2.stdout,stdout=PIPE)
  p4 = Popen(cmd4.split(),stdin=p3.stdout,stdout=PIPE)
  p5 = Popen(cmd5.split(),stdin=p4.stdout,stdout=PIPE)
  res = p5.communicate()[0].rstrip()
  p3.communicate()
  p4.communicate()
  p2.communicate()
  p1.communicate()
  inf.close()
  #cmd = "samtools view "+fname+" | cut -f 1 | sort | uniq | wc -l"
  #p = Popen(cmd,shell=True,stdout=PIPE)
  #res =p.communicate()[0]
  return [fname,res]

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inputs',nargs='+',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  args = parser.parse_args()
  # Temporary working directory step 2 of 3 - Creation
  return args

if __name__=="__main__":
  main()
