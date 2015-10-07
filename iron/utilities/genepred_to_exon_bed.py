#!/usr/bin/python

import sys, re, os, argparse
from GenePredBasics import line_to_entry as genepred_line_to_entry


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help='FILENAME input genepred, use - for STDIN')
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  for line in inf:
    if re.match('^#',line): continue
    e = genepred_line_to_entry(line)
    for i in range(0,len(e['exonStarts'])):
      print e['chrom']+"\t"+str(e['exonStarts'][i])+"\t"+str(e['exonEnds'][i])+"\t"+e['gene_name']+"\t"+e['name']+"\t"+str(i)
  inf.close()


if __name__=="__main__":
  main()
