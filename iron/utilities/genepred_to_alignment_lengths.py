#!/usr/bin/python
import argparse, sys
import GenePredBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="GENEPRED file input use - for STDIN")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  for line in inf:
    e = GenePredBasics.GenePredEntry()
    e.line_to_entry(line.rstrip())
    print e.entry['gene_name']+"\t"+e.entry['name']+"\t"+str(e.length())
  inf.close()

if __name__=="__main__":
  main()
