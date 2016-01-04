#!/usr/bin/python
import argparse, sys
from PSLBasics import PSL


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="PSLFILE or - for STDIN")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  z = 0
  for line in inf:
    z+=1
    p = PSL(line.rstrip())
    print str(z) + "\t" + p.value('qName') + "\t" + p.value('tName')+"\t"+str(p.get_coverage())+"\t"+str(p.value('qSize'))+"\t"+str(p.get_quality())
  inf.close()

if __name__=="__main__":
  main()
