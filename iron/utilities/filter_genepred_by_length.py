#!/usr/bin/python
import argparse, sys, re
from GenePredBasics import line_to_entry

def main():
  parser = argparse.ArgumentParser(description="Filter a genepred by transcript length")
  parser.add_argument('-i','--input',help="Input '-' for STDOUT")
  parser.add_argument('--min',type=int,help="Minimum transcript length")
  parser.add_argument('--max',type=int,help="Maximum transcript length")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  for line in inf:
    if re.match('^#',line): continue
    line = line.rstrip()
    e = line_to_entry(line)
    tot = 0
    for i in range(0,len(e['exonStarts'])):
      tot += e['exonEnds'][i]-e['exonStarts'][i]
    if args.min:
      if tot < args.min:
        continue
    if args.max:
      if tot > args.max:
        continue
    # If we are still here we can print
    print line
  
if __name__=="__main__":
  main()
