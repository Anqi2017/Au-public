#!/usr/bin/python
import argparse, sys, re
import PSLBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('psl_file',help="use - for STDIN")
  parser.add_argument('--filter',action='store_true',help="Only output passing lines")
  args = parser.parse_args()
  inf = sys.stdin
  if args.psl_file != '-':
    inf = open(args.psl_file)
  z = 0
  for line in inf:
    z += 1
    if args.filter:
      if PSLBasics.is_valid(line): print line.rstrip()
      continue
    if not PSLBasics.is_valid(line):
      print "bad line "+str(z)
      print line
      return
  if not args.filter:
    print "PSL file looks good"
  inf.close()

if __name__=="__main__":
  main()
