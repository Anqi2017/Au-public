#!/usr/bin/python
import sys, argparse
from Bio.Format.GPD import GPDStream
from Bio.Stream import LocusStream
from Bio.Range import ranges_to_coverage

def main():
  parser = argparse.ArgumentParser(description="Convert sorted gpd file to bed depth",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  args = parser.parse_args()
  
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)
  loci = LocusStream(GPDStream(args.input))
  for locus in loci:
    exranges = []
    for entry in locus.get_payload():
      for exon in entry.exons:
        exranges.append(exon.get_range())
    covs = ranges_to_coverage(exranges)    
    for cov in covs:    
      print "\t".join([str(x) for x in cov.get_bed_coordinates()])+"\t"+str(+cov.get_payload())
if __name__=="__main__":
  main()
