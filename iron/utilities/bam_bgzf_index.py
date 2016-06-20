#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree, copy
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from Bio.Format.Sam import BAMFile

def main(args):
  ind_path = args.input+'.bgi'
  if args.output: ind_path = args.output
  if os.path.isfile(ind_path) and not args.output:
    sys.stderr.write("ERROR index file already there.  Delete it if you want to rebuild a new one.\n")
    sys.exit()
  bf = BAMFile(args.input)
  bf.write_index(args.tempdir+'/myfile.bgi',verbose=True)
  copy(args.tempdir+'/myfile.bgi',ind_path)  

  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Generate our .bgi index for a bam file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT BAM FILE")
  parser.add_argument('--output','-o',help="Specifiy path to write index")
  #parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  group.add_argument('--no_primary_search',action='store_true',help="Dont try to assign a primary flag")
  args = parser.parse_args()
  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

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

def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main(args)
