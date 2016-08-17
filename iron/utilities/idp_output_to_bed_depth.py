#!/usr/bin/python
import sys, argparse, gzip
from Bio.Format.GPD import GPD
from Bio.Range import ranges_to_coverage

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="IDP output folder")
  #parser.add_argument('--min_exons',type=int,default=1,help="At least this number of exons")
  parser.add_argument('--offset',type=int,default=1,help="add this much to all expressions")
  parser.add_argument('--mult',type=int,default=10,help="multiply all expressions by this much")
  parser.add_argument('-o','--output',help="OUTPUT file or nothing for STDOUT")
  args = parser.parse_args()
  
  args.input= args.input.rstrip('/')
  inf = open(args.input+'/isoform.gpd')
  sys.stderr.write("Reading isoform.gpd\n")
  txs = {}
  for line in inf:
    gpd = GPD(line)
    tx = gpd.get_transcript_name()
    if tx not in txs:
      txs[tx] = []
    for exon in gpd.exons:
      txs[tx].append(exon.get_range())
  inf.close()

  sys.stderr.write("Reading isoform.exp file\n")
  inf = open(args.input+'/isoform.exp')
  vals = []
  for line in inf:
      f = line.rstrip().split("\t")
      v = int((float(f[1])*args.mult)+args.offset)
      tx = f[0]
      exons = txs[tx]
      #if len(exons) < args.min_exons: continue
      for i in range(0,v):
        vals += exons[:]
  inf.close()
  sys.stderr.write("Generating coverage file "+str(len(vals))+"\n")
  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz':
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  covs = ranges_to_coverage(vals)
  for v in covs:
    of.write(v.chr+"\t"+str(v.start-1)+"\t"+str(v.end)+"\t"+str(v.get_payload())+"\n")
  #    of.write(tx+"\t"+gene+"\t"+str(genes[gene]['transcripts'][tx])+"\t"+str(genes[gene]['cnt'])+"\n")
  of.close()
if __name__=="__main__":
  main()
