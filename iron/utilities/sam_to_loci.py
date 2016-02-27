#!/usr/bin/python
import sys, argparse
from SamBasics import SamLocusStream

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  args = parser.parse_args()
  
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)
  reader = SamLocusStream(args.input)
  while True:
    locus = reader.read_locus()
    if not locus: break
    print locus[0].chr+"\t"+str(locus[0].start-1)+"\t"+str(locus[0].end)
if __name__=="__main__":
  main()
