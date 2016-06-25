#!/usr/bin/python
import argparse, os, inspect, sys
from subprocess import Popen, PIPE
from tempfile import mkdtemp, gettempdir
from multiprocessing import cpu_count
from shutil import rmtree

#bring in the folder to the path for our utilities
pythonfolder_loc = "../utilities"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import prepare_all_data
import create_html

g_version = None

def main(args):

  if not args.output and not args.portable_output and not args.output_folder:
    sys.stderr.write("ERROR: must specify some kind of output\n")
    sys.exit()

  ## Check and see if directory for outputs exists
  if args.output_folder:
    if os.path.isdir(args.output):
      sys.stderr.write("ERROR: output directory already exists.  Remove it to to use this location\n")
      sys.exit()


  global g_version
  #Make sure rscript is installed
  try:
    cmd = 'Rscript --version'
    prscript = Popen(cmd.split(),stdout=PIPE,stderr=PIPE)
    rline = prscript.communicate()
    sys.stderr.write("Using Rscript version:\n")
    sys.stderr.write(rline[1].rstrip()+"\n")
  except:
    sys.stderr.write("ERROR: Rscript not installed\n")
    sys.exit()

  # Check and see if directory for output exists
  if not args.output and not args.portable_output:
    sys.stderr.write("ERROR: Must specify some kind of output\n")
    sys.exit()

  if args.no_reference:
    sys.stderr.write("WARNING: No reference specified.  Will be unable to report error profile\n")
  if args.no_annotation:
    sys.stderr.write("WARNING: No annotation specified.  Will be unable to report feature specific outputs\n")

  prepare_all_data.external(args)
  create_html.external(args,version=g_version)

  # Temporary working directory step 3 of 3 - Cleanup
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

def external_cmd(cmd,version=None):
  #set version by input
  global g_version
  g_version = version

  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

def do_inputs():
  parser=argparse.ArgumentParser(description="Create an output report",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")

  # output options
  parser.add_argument('-o','--output',help="OUTPUT xhtml with data")
  parser.add_argument('--portable_output',help="OUTPUT file in a small xhtml format")
  parser.add_argument('--output_folder',help="OUTPUT folder of all data")

  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('-r','--reference',help="Reference Fasta")
  group1.add_argument('--no_reference',action='store_true',help="No Reference Fasta")
  group2 = parser.add_mutually_exclusive_group(required=True)
  group2.add_argument('--annotation',help="Reference annotation genePred")
  group2.add_argument('--no_annotation',action='store_true',help="No annotation is available")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")

  ### Parameters for alignment plots
  parser.add_argument('--min_aligned_bases',type=int,default=50,help="for analysizing alignment, minimum bases to consider")
  parser.add_argument('--max_query_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  parser.add_argument('--max_target_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  parser.add_argument('--max_query_gap',type=int,help="for testing gapped alignment advantge")
  parser.add_argument('--max_target_gap',type=int,default=500000,help="for testing gapped alignment advantage")
  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="require gapped alignment to be this much better (in alignment length) than single alignment to consider it.")
  
  ### Parameters for locus analysis
  parser.add_argument('--do_loci',action='store_true',help="this analysis is time consuming at the moment\n")
  parser.add_argument('--min_depth',type=float,default=1.5,help="require this or more read depth to consider locus")
  parser.add_argument('--min_coverage_at_depth',type=float,default=0.8,help="require at leas this much of the read be covered at min_depth")
  parser.add_argument('--min_exon_count',type=int,default=2,help="Require at least this many exons in a read to consider assignment to a locus")
  parser.add_argument('--locus_downsample',type=int,default=100,help="Limit how deep to search loci\n")

  ### Params for alignment error plot
  parser.add_argument('--alignment_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  parser.add_argument('--alignment_error_max_length',type=int,default=100000,help="The maximum number of alignment bases to calculate error from")
  
  ### Params for context error plot
  parser.add_argument('--context_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  parser.add_argument('--context_error_stopping_point',type=int,default=1000,help="Sample at least this number of each context")

  ## Params for rarefraction plots
  parser.add_argument('--samples_per_xval',type=int,default=500)

  args = parser.parse_args()
  setup_tempdir(args)
  return args

if __name__=='__main__':
  do_inputs()
  main(args)
