#!/usr/bin/python
import argparse,sys
from PSLBasics import PSL

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="Use - for STDIN")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--psl',action='store_true')
  group.add_argument('--gpd',action='store_true')
  group.add_argument('--fasta',action='store_true')
  group.add_argument('--fastq',action='store_true')
  args = parser.parse_args()
  
  if args.input == '-':  args.input = sys.stdin
  else: args.input = open(args.input)

  if args.psl:  do_psl(args)

def do_psl(args):
  for line in args.input:
    psl = PSL(line)
    cov = sum(psl.value('blockSizes'))
    print cov

if __name__=="__main__":
  main()
