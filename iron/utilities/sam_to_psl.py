#!/usr/bin/python
import sys, argparse, inspect, os

#bring in the folder to the path for our utilities
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

from Bio.Format.Fasta import FastaData
from Bio.Format.Sam import BAMFile, SamStream

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN or specify a BAM file")
  parser.add_argument('-r','--reference',help="Reference fasta",required=True)
  args = parser.parse_args()

  ref = None
  if args.reference:
    ref = FastaData(open(args.reference,'rb').read())
  
  if args.input == '-':
    args.input = SamStream(sys.stdin,reference=ref)
  else: args.input = BAMFile(args.input,reference=ref)
  for e in args.input:
    if e.is_aligned():
      print e.get_PSL()

if __name__=="__main__":
  main()
