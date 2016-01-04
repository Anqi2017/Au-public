#!/usr/bin/python
import sys, re, argparse
from SequenceBasics import read_fasta_into_hash

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('ref_genome')
  parser.add_argument('phased_vcf')
  args = parser.parse_args()
  g = read_fasta_into_hash(args.ref_genome)
  gL = {}
  gR = {}
  for chr in g:
    gL[chr] = [x for x in g[chr].upper()]
    gR[chr] = [x for x in g[chr].upper()]
  #with open('1KG.biallelic.het.exonic/1KG.biallelic.het.exonic.vcf') as inf:
  z = 0
  with open(args.phased_vcf) as inf:
    for line in inf:
      if re.match('#',line): continue
      z += 1
      sys.stderr.write(str(z)+"\r")
      f = line.rstrip().split("\t")
      chr = f[0]
      [n1,n2] = [int(x) for x in f[9].split('|')]
      if int(f[1]) > len(gL[chr]): 
        sys.stderr.write(line)
        sys.exit()
      if n1 == 0: 
        gL[chr][int(f[1])-1] = f[3]
      else: 
        gL[chr][int(f[1])-1] = f[4]
      if n2 == 0:
        gR[chr][int(f[1])-1] = f[3]
      else:
        gR[chr][int(f[1])-1] = f[4]
  sys.stderr.write("\nalmost done\n")
  ofL = open('L.fa','w')
  for chr in sorted(gL.keys()):
    ofL.write(">"+chr+"\n")
    ofL.write("".join(gL[chr])+"\n")
  ofL.close()
  ofR = open('R.fa','w')
  for chr in sorted(gR.keys()):
    ofR.write(">"+chr+"\n")
    ofR.write("".join(gR[chr])+"\n")
  ofR.close()
if __name__=="__main__":
  main()
