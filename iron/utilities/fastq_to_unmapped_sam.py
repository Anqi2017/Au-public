#!/usr/bin/python
import sys, argparse
from Bio.Format.Fasta import FastaData
from Bio.Format.Fastq import FastqHandle

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('-r','--reference',help="reference genome FASTA")
  parser.add_argument('--no_qual',action='store_true',help="dont put in quality")
  
  args = parser.parse_args()
  ref = {}
  if args.reference:
    ref = FastaData(open(args.reference,'rb').read())
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)
  
  
  h1 = '@HD	VN:1.0	SO:unsorted'
  h2 = '@PG	ID:FA2UN	PN:FA2UN	VN:2016-06-09	CL:'+' '.join(sys.argv)
  print h1
  print h2
  if ref:
    for chr in sorted(ref.keys()):
      print "@SQ\tSN:"+chr+"\t"+'LN:'+str(len(ref[chr]))
  inf = FastqHandle(args.input)
  for e in inf:
    o =  ''
    o += e.name+"\t"
    o += "4\t"
    o += "*\t"
    o += "0\t"
    o += "0\t"
    o += "*\t"
    o += "*\t"
    o += "0\t"
    o += "0\t"
    o += e.seq+"\t"
    if args.no_qual:
      o+= "*\t"
    else:
      o += e.qual+"\t"
    o += "XO:Z:NM"
    print o

if __name__=="__main__":
  main()
