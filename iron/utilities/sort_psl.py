#!/usr/bin/python
import argparse, sys, os, random, subprocess
from shutil import rmtree

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("psl_file",help="PSLFILE or - for sam streamed to stdin")
  parser.add_argument("-o","--output",help="OUTFILE or STDOUT if not specified")
  parser.add_argument('-T',"--tempdir",help="DIR to store a temp directory")
  parser.add_argument('-S',"--size",help="S option for linux sort units are kb unless specified")
  args = parser.parse_args()

  # Manage your tempdir
  # put it into args.tempdir
  if args.tempdir:
    rnum = random.randint(1,100000000)
    args.tempdir = args.tempdir.rstrip('/')+'/weirathe.'+str(rnum) 
    if not os.path.exists(args.tempdir):
      os.makedirs(args.tempdir)
  
  inf = sys.stdin
  if args.psl_file != '-':
    inf = open(args.psl_file)
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')

  cmd = "sort -k 10,10"
  if args.tempdir:
    cmd += " -T "+args.tempdir
  if args.size:
    cmd += " -S "+args.size
  p1 = subprocess.Popen(cmd.split(),stdin = subprocess.PIPE, stdout = of)
  for line in inf:
    p1.stdin.write(line)
  p1.communicate()
  inf.close()
  of.close()
  #cleanup
  if args.tempdir:
    rmtree(args.tempdir)

if __name__=="__main__":
  main()
