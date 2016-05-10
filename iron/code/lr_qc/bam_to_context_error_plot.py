#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from Bio.Format.Sam import BAMFile
from Bio.Errors import ErrorProfileFactory
from Bio.Format.Fasta import FastaData
from subprocess import call


def main():
  #do our inputs
  args = do_inputs()
  # make our error profile report
  sys.stderr.write("Reading reference fasta\n")
  ref = FastaData(open(args.reference).read())
  sys.stderr.write("Reading alignments\n")
  bf = BAMFile(args.input,reference=ref)
  epf = ErrorProfileFactory()
  z = 0
  strand = 'target'
  if args.query: strand = 'query'
  con = 0
  for e in bf:
    if e.is_aligned():
      epf.add_alignment(e)
      z+=1
      if z%100==1:
        con = epf.get_min_context_count(strand)
      sys.stderr.write(str(z)+" alignments, "+str(con)+" min context coverage\r")
      if args.max_alignments <= z: break
      if args.stopping_point <= con: break
  sys.stderr.write("\n")
  print 'working with:'
  print str(z)+" alignments, "+str(con)+" min context coverage"
  epf.write_context_error_report(args.tempdir+'/err.txt',strand)  

  cmd = os.path.dirname(os.path.realpath(__file__))+'/plot_base_error_context.r '+args.tempdir+'/err.txt '+args.output
  print cmd
  call(cmd.split())
  print 'finished'
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-r','--reference',required=True,help="Reference Genome")
  parser.add_argument('-o','--output',required=True,help="OUTPUTFILE")

  parser.add_argument('--max_alignments',type=int,default=10000000000,help="The maximum number of alignments to scan")
  parser.add_argument('--stopping_point',type=int,default=1000,help="Stop after you see this many of each context")
  #parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('--target',action='store_true',help="Context on the target sequence")
  group1.add_argument('--query',action='store_true',help="Context on the query sequence")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
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

if __name__=="__main__":
  main()
