#!/usr/bin/python
import argparse, sys, re
from GenePredBasics import GenePredEntry as GPD

def main():
  parser = argparse.ArgumentParser(description="Filter a genepred by transcript length")
  parser.add_argument('input',help="Input '-' for STDOUT")
  parser.add_argument('--min_length',type=int,help="Minimum transcript length")
  parser.add_argument('--max_length',type=int,help="Maximum transcript length")
  parser.add_argument('--names',help="filter on a name list")
  parser.add_argument('--gene_names',help="filter on a gene name list")
  parser.add_argument('-v','--invert',action='store_true',help='Invert search result')
  args = parser.parse_args()
  name_list = set()
  gene_name_list = set()
  if args.names:
    with open(args.names) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        name_list.add(f[0])
  if args.gene_names:
    with open(args.gene_names) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        gene_name_list.add(f[0])
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  for line in inf:
    if re.match('^#',line): continue
    is_good = True
    g = GPD(line.rstrip())
    tot = g.length()
    if args.min_length:
      if tot < args.min_length:
        is_good = False
    if args.max_length:
      if tot > args.max_length:
        is_good = False
    if args.names:
      if g.value('name') not in name_list:
        is_good = False
    if args.gene_names:
      if g.value('gene_name') not in args.gene_name_list:
        is_good = False
    # If we are still here we can print
    if not args.invert:
      if is_good: print line.rstrip()
    else:
      if not is_good: print line.rstrip()
  
if __name__=="__main__":
  main()
