#!/usr/bin/python
import argparse, sys, os, random
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from Bio.Format.Sam import BAMFile
from Bio.Errors import ErrorProfileFactory
from Bio.Format.Fasta import FastaData
from subprocess import call

# Take the bam file as an input and produce plots and data file for context errors.

def main():
  #do our inputs
  args = do_inputs()
  # make our error profile report
  sys.stderr.write("Reading reference fasta\n")
  ref = FastaData(open(args.reference).read())
  sys.stderr.write("Reading alignments\n")
  epf = ErrorProfileFactory()
  if args.random:
    bf = BAMFile(args.input,reference=ref)
    bf.read_index()
    if not bf.has_index():
      sys.stderr.write("Random access requires an index be set\n")
    z = 0
    strand = 'target'
    if args.query: strand = 'query'
    con = 0
    while True:
      rname = random.choice(bf.index.get_names())
      #print rname
      coord = bf.index.get_longest_target_alignment_coords_by_name(rname)
      #print coord
      e = bf.fetch_by_coord(coord)
      if e.is_aligned():
        epf.add_alignment(e)
        z+=1
        if z%100==1:
          con = epf.get_min_context_count(strand)
        sys.stderr.write(str(z)+" alignments, "+str(con)+" min context coverage\r")
        if args.max_alignments <= z: break
        if args.stopping_point <= con: break
    
  else:
    bf = BAMFile(args.input,reference=ref)
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
  sys.stderr.write('working with:'+"\n")
  sys.stderr.write(str(z)+" alignments, "+str(con)+" min context coverage"+"\n")
  epf.write_context_error_report(args.tempdir+'/err.txt',strand)  

  for ofile in args.output:
    cmd = 'Rscript '+os.path.dirname(os.path.realpath(__file__))+'/plot_base_error_context.r '+args.tempdir+'/err.txt '+ofile+' '
    if args.scale:
      cmd += ' '.join([str(x) for x in args.scale])
    sys.stderr.write(cmd+"\n")
    call(cmd.split())
  sys.stderr.write("finished\n")
  if args.output_raw:
    of = open(args.output_raw,'w')
    with open(args.tempdir+"/err.txt") as inf:
      for line in inf:
        of.write(line)
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-r','--reference',required=True,help="Reference Genome")
  parser.add_argument('-o','--output',nargs='+',required=True,help="OUTPUTFILE(s)")
  parser.add_argument('--output_raw',help="Save the raw data")
  parser.add_argument('--scale',type=float,nargs=6,help="<insertion_min> <insertion_max> <mismatch_min> <mismatch_max> <deletion_min> <deletion_max>")
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
  parser.add_argument('--random',action='store_true',help="Randomly select alignments, requires an indexed bam")
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
