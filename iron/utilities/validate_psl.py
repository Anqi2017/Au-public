#!/usr/bin/python
import argparse, sys, re
import PSLBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('psl_file',help="use - for STDIN")
  args = parser.parse_args()
  inf = sys.stdin
  if args.psl_file != '-':
    inf = open(args.psl_file)
  z = 0
  for line in inf:
    z += 1
    if not PSLBasics.is_valid(line):
      print "bad line "+str(z)
      print line
      return
  print "PSL file looks good"
  inf.close()

if __name__=="__main__":
  main()
