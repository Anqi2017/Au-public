#!/usr/bin/python
import sys, argparse, re, gzip
from Bio.Format.GPD import GPDStream

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="best gpd")
  parser.add_argument('-o','--output',help="write to output")
  args = parser.parse_args()
  
  inf = None
  if re.search('\.gz$',args.input):
    inf = gzip.open(args.input)
  else:
    inf = open(args.input)
  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')

  gs = GPDStream(inf)
  for gpd in gs:
    of.write(str(gpd.get_length())+"\t"+str(gpd.get_exon_count())+"\n")
  of.close()  


if __name__=="__main__":
  main()
