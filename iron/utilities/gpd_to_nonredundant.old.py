#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree
from multiprocessing import cpu_count, Pool, Lock
import GenePredBasics, GenePredFuzzyBasics

glock = Lock()
gout = sys.stdout

def main():
  #do our inputs
  args = do_inputs()
  global gout
  gout = args.output
  gls = GenePredBasics.GenePredLocusStream(args.input)
  fgs = GenePredFuzzyBasics.FuzzyGenePredSeparator()
  if args.threads > 1:
    p = Pool(processes=args.threads)
  while True:
    buffer = gls.read_locus()
    if not buffer: break
    if args.threads > 1:
      p.apply_async(process_buffer,args=(buffer,args),callback=out_gpds)
    else:
      v = process_buffer(buffer,args)
      out_gpds(v)
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\n")

def process_buffer(buffer,args):
  fzs = GenePredFuzzyBasics.greedy_gpd_list_to_combined_fuzzy_list(buffer,args.junction_tolerance)
  return [fzs,args]

def out_gpds(v):
  global glock
  global gout
  glock.acquire()
  [fzs,args] = v
  line = None
  for fz in fzs:
    line = fz.get_genepred_line()
    gout.write(line+"\n")
  if line:
    f = line.rstrip().split("\t")
    status = f[2]+":"+str(int(f[4])+1)+'-'+f[5]
    sys.stderr.write("\r"+status+"         ")
  glock.release()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="")
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  parser.add_argument('--junction_tolerance',type=int,default=10)
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
  #setup_tempdir(args)
  return args

if __name__=="__main__":
  main()
