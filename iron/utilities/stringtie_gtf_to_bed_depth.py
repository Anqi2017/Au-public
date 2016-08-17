#!/usr/bin/python
import sys, argparse, gzip, re
from Bio.Format.GPD import GPD
from Bio.Range import ranges_to_coverage, Bed

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--offset',type=int,default=0,help="add this much to transcript tpms")
  parser.add_argument('--mult',type=int,default=10,help="multiply this much to tpms")
  parser.add_argument('--min_exons',type=int,default=1,help="require at least this many exons")
  parser.add_argument('-o','--output',help="OUTPUT file or nothing for STDOUT")
  args = parser.parse_args()
  
  inf = sys.stdin
  if args.input != '-':
    if args.input[-3:]=='.gz':
      inf = gzip.open(args.input)
    else: inf = open(args.input)
  genes = {}
  sys.stderr.write("Reading gtf file\n")
  txs = {}
  for line in inf:
    if re.match('#',line): continue
    f = line.rstrip().split("\t")
    tx = None
    if f[2] == 'exon' or f[2] == 'transcript':
      tx = re.search('transcript_id\s+"([^"]+)"',f[8]).group(1)
      if tx not in txs:
        txs[tx] = {}
        txs[tx]['tpm'] = 0
        txs[tx]['exons'] = []
    if f[2] == 'transcript':
      tpm = float(re.search('TPM\s+"([^"]+)"',f[8]).group(1))
      txs[tx]['tpm'] = int((tpm*float(args.mult))+args.offset)
    if f[2] == 'exon':
      chr = f[0]
      start = int(f[3])-1
      end = int(f[4])
      txs[tx]['exons'].append(Bed(chr,start,end))
  inf.close()
  vals = []
  sys.stderr.write("Traversing annotation file\n")
  for tx in txs:
    exons = txs[tx]['exons']
    v = txs[tx]['tpm']
    if len(exons) < args.min_exons: continue
    for i in range(0,v):
      vals += exons[:]
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
