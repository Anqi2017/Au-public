#!/usr/bin/python
import sys, argparse, re, gzip
from Bio.Format.GPD import GPDStream
from Bio.Stream import LocusStream
from Bio.Range import ranges_to_coverage, union_range_array

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('-o','--output',help="output file or use STDOUT if not set")
  args = parser.parse_args()
  
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)
  gs = GPDStream(args.input)
  ls = LocusStream(gs)
  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output): 
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  for rng in ls:
    sys.stderr.write(rng.get_range_string()+"    \r")
    gpds = rng.get_payload()
    exs = []
    for ex_set in [[y.get_range() for y in x.exons] for x in gpds]:
      exs += ex_set
    cov = ranges_to_coverage(exs)
    #use our coverage data on each gpd entry now
    for gpd in gpds:
      totcov = 0
      for exon in [x.get_range() for x in gpd.exons]:
        gcovs = union_range_array(exon,cov,payload=2)
        totcov += sum([x.get_payload()*x.length() for x in gcovs])
      of.write(gpd.get_gene_name()+"\t"+str(gpd.get_exon_count())+"\t"+str(gpd.get_length())+"\t"+str(float(totcov)/float(gpd.get_length()))+"\n")
  sys.stderr.write("\n")
  of.close()
if __name__=="__main__":
  main()
