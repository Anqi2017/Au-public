#!/usr/bin/python
import argparse, sys
import GenePredBasics

def main():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('genepred',help="FILENAME or use - for STDIN")
  parser.add_argument('--smoothing_size',type=int,default=68,help="INT no gaps less than this size")
  args = parser.parse_args()
  inf = sys.stdin
  if args.genepred != '-':
    inf = open(args.genepred)
  for line in inf:
    e = GenePredBasics.line_to_entry(line)
    e2 = GenePredBasics.smooth_gaps(e,args.smoothing_size)
    print GenePredBasics.entry_to_line(e2)

if __name__=="__main__":
  main()
