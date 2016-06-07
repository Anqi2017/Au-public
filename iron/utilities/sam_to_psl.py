#!/usr/bin/python
import sys, argparse
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
