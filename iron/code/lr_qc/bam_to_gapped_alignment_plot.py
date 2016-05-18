#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree, copy
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import call

def main():
  #do our inputs
  args = do_inputs()
  udir = os.path.dirname(os.path.realpath(__file__))
  sys.stderr.write("Making text report\n")

  cmd = 'python '+udir+'/bam_to_gapped_alignment_report.py '+args.input+' -o '+args.tempdir+'/report.txt'+' '
  if args.max_query_overlap:
    cmd += '--max_query_overlap '+str(args.max_query_overlap)+' '
  if args.max_target_overlap:
    cmd += '--max_target_overlap '+str(args.max_target_overlap)+' '
  if args.max_target_gap:
    cmd += '--max_target_gap '+str(args.max_target_gap)+' '
  if args.max_query_gap:
    cmd += '--max_query_gap '+str(args.max_query_gap)+' '
  if args.required_fractional_improvement:
    cmd += '--required_fractional_improvement '+str(args.required_fractional_improvement)+' '
  sys.stderr.write(cmd+"\n")
  call(cmd.split())

  sys.stderr.write("Finished making report\n")
  sys.stderr.write("making plot\n")
  for ofile in args.output:
    cmd = 'Rscript '+udir +'/plot_gapped_alignment_statistics.r '+args.tempdir+'/report.txt '+ofile
    sys.stderr.write(cmd+"\n")
    call(cmd.split())

  if args.output_raw:
    copy(args.tempdir+"/report.txt",args.output_raw)
  if args.output_stats:
    do_stats(args)
  sys.stderr.write("Finished.\n")
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_stats(args):
  total_reads = 0
  unaligned_reads = 0
  aligned_reads = 0
  single_align_reads = 0
  gapped_align_reads = 0
  total_bases = 0
  unaligned_bases = 0
  aligned_bases = 0
  single_align_bases = 0
  gapped_align_bases = 0
  with open(args.tempdir+"/report.txt") as inf:
    for line in inf:
      (type, single, both, rlen) = [int(x) for x in line.rstrip().split("\t")]
      total_reads += 1
      if type==0: unaligned_reads += 1
      else: aligned_reads += 1
      if type==1: single_align_reads +=1
      if type>1: gapped_align_reads += 1
      if type>0:
        total_bases += rlen
        unaligned_bases += (rlen-both)
        aligned_bases += both
        single_align_bases += single
        gapped_align_bases += both-single
  of = open(args.output_stats,'w')
  of.write("TOTAL_READS\t"+str(total_reads)+"\n")
  of.write("UNALIGNED_READS\t"+str(unaligned_reads)+"\n")
  of.write("ALIGNED_READS\t"+str(aligned_reads)+"\n")
  of.write("SINGLE_ALIGN_READS\t"+str(single_align_reads)+"\n")
  of.write("GAPPED_ALIGN_READS\t"+str(gapped_align_reads)+"\n")
  of.write("TOTAL_BASES\t"+str(total_bases)+"\n")
  of.write("UNALIGNED_BASES\t"+str(unaligned_bases)+"\n")
  of.write("ALIGNED_BASES\t"+str(aligned_bases)+"\n")
  of.write("SINGLE_ALIGN_BASES\t"+str(single_align_bases)+"\n")
  of.write("GAPPED_ALIGN_BASES\t"+str(gapped_align_bases)+"\n")
  of.close()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT BAMFILE")
  parser.add_argument('-o','--output',nargs='+',help="OUTPUT FILE can put multiple")
  parser.add_argument('--output_raw',help="Save the raw data here")
  parser.add_argument('--output_stats',help="Save some summary statistics")
  parser.add_argument('--reference','-r',required=True,help="Fasta reference file")

  # Args for the gapped alignment report
  parser.add_argument('--max_query_overlap',type=int,default=10,help="Consider two alignments incompatible if greater than this")
  parser.add_argument('--max_target_overlap',type=int,default=10,help="Consider two alignments incompatible if greater than this")
  parser.add_argument('--max_target_gap',type=int,default=500000,help="Consider a gapped alignment incompatible if greater than this")
  parser.add_argument('--max_query_gap',type=int,help="Consider a gapped alignment incompatible if greater thant this")
  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="Result should be this much better than the original")

  # Temporary working directory step 1 of 3 - Definition
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
