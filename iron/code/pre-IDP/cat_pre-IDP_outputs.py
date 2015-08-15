#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree, copyfile, copytree

def main():
  parser = argparse.ArgumentParser(description="Take outputs from any of the pre-IDP ",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--stage',type=int,choices=[1,2,3],required=True,help='INT stage')
  parser.add_argument('--output',required=True,help='FOLDERNAME for output')
  parser.add_argument('folder_names',nargs='+',help='FOLDERNAME(s)')
  args = parser.parse_args()
  if os.path.exists(args.output):
    sys.stderr.write("ERROR: output directory already exists\n")
    sys.exit()
  args.output=args.output.rstrip('/')
  os.makedirs(args.output)
  if args.stage == 1:
    process_stage_1(args)

def process_stage_1(args):
  os.makedirs(args.output+'/ccs_hq')
  os.makedirs(args.output+'/ccs_lq')
  os.makedirs(args.output+'/subreads')
  of1 = open(args.output+'/ccs_lq_and_longest_subreads_to_correct.fa','w')
  of2 = open(args.output+'/lr_nonredundant_uncorrected.fa','w')
  of3 = open(args.output+'/LOG','w')
  of4 = open(args.output+'/ccs_hq/ccs_hq.fa','w')
  of5 = open(args.output+'/ccs_hq/ccs_hq.fq','w')
  of6 = open(args.output+'/ccs_lq/ccs_lq.fa','w')
  of7 = open(args.output+'/ccs_lq/ccs_lq.fq','w')
  of8 = open(args.output+'/subreads/subreads.fa','w')
  of9 = open(args.output+'/subreads/subreads.fq','w')
  for folder in args.folder_names:
    if not os.path.isdir(folder):
      sys.stderr.write("ERROR: given folder name "+folder+" is not a folder\n")
      sys.exit
    with open(folder.rstrip('/')+'/ccs_lq_and_longest_subreads_to_correct.fa') as inf:
      for line in inf:
        of1.write(line)
    with open(folder.rstrip('/')+'/lr_nonredundant_uncorrected.fa') as inf:
      for line in inf:
        of2.write(line)
    with open(folder.rstrip('/')+'/LOG') as inf:
      for line in inf:
        of3.write(line)
    with open(folder.rstrip('/')+'/ccs_hq/ccs_hq.fa') as inf:
      for line in inf:
        of4.write(line)
    with open(folder.rstrip('/')+'/ccs_hq/ccs_hq.fq') as inf:
      for line in inf:
        of5.write(line)
    with open(folder.rstrip('/')+'/ccs_lq/ccs_lq.fa') as inf:
      for line in inf:
        of6.write(line)
    with open(folder.rstrip('/')+'/ccs_lq/ccs_lq.fq') as inf:
      for line in inf:
        of7.write(line)
    with open(folder.rstrip('/')+'/subreads/subreads.fa') as inf:
      for line in inf:
        of8.write(line)
    with open(folder.rstrip('/')+'/subreads/subreads.fq') as inf:
      for line in inf:
        of9.write(line)
  of1.close()
  of2.close()
  of3.close()
  of4.close()
  of5.close()
  of6.close()
  of7.close()
  of8.close()
  of9.close()
  return

if __name__ == "__main__":
  main()
