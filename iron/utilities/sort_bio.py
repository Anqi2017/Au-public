#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree
#from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="")
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--name',action='store_true',help="Sort by query name rather than location.  For GenePred this will default to gene name then the transcript name.")
  #parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--gpd',action='store_true')
  group2.add_argument('--bed',action='store_true')
  group2.add_argument('--psl',action='store_true')
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  # Setup inputs 
  if args.input == '-':
    args.input = sys.stdin
  else:
    args.input = open(args.input)
  # Setup outputs
  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout
  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

def main():
  #do our inputs
  args = do_inputs()
  # Temporary working directory step 3 of 3 - Cleanup
  #sys.stderr.write("working in: "+args.tempdir+"\n")
  cmd = "sort -T "+args.tempdir+'/'
  if args.psl:
    if args.name:
      cmd = "sort -k10,10 -T "+args.tempdir+'/'
    else:
      cmd = "sort -k14,14 -k15,15n -k16,16n -k9,9 -T "+args.tempdir+'/'
  if args.bed:
    cmd = "sort -k1,1 -k2,2n -T "+args.tempdir+'/'
  if args.gpd:
    if args.name:
      cmd = "sort -k1,1 -k2,2 -T "+args.tempdir+'/'
    else:
      cmd = "sort -k3,3 -k5,5n -k6,6n -k4,4 -T "+args.tempdir+'/'
  p = Popen(cmd.split(),stdout=args.output,stdin=PIPE)
  for line in args.input:
    p.stdin.write(line)
  p.stdin.close()
  p.wait()
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

if __name__=="__main__":
  main()
