#!/usr/bin/python
import sys, argparse
from subprocess import PIPE, Popen
from multiprocessing import cpu_count, Pool, Lock

glock = Lock()
of = None
def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',nargs='+',help="Files to do")
  parser.add_argument('-o','--output',help="output or - for stdout")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="number of processes")
  args = parser.parse_args()
  global of
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for file in args.input:
    if args.threads > 1:
      p.apply_async(do_file,args=(file,),callback=do_out)
    else:
      r = do_file(file)
      do_out(r)
  if args.threads > 1:
    p.close()
    p.join()
  of.close()
def do_out(r):
  global glock
  glock.acquire()
  global of
  of.write(r+"\n")
  glock.release()

def do_file(file):
  cmd = "du -hs "+file
  p = Popen(cmd.split(),stdout=PIPE)
  return p.communicate()[0].rstrip()
  

if __name__=="__main__":
  main()
