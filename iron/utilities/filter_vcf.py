#!/usr/bin/python
import argparse, sys, re
import VCFBasics
import re

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input_vcf',help="VCFFILE or - for STDIN")
  parser.add_argument('--quality',type=float,help="FLOAT minimum quality")
  #parser.add_argument('--biallelic',action='store_true')
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--mintotaldepth',type=int,help="Require between all alleles")
  group.add_argument('--mindepth',type=int,help="Require this between the two strands for the lower expressed allele")
  group.add_argument('--mindepthstrand',type=int,help="Require this depth on both strands")
  parser.add_argument('--biallelic',action='store_true',help="Require line type be biallelic")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input_vcf != '-': inf = open(args.input_vcf)
  for line in inf:
    v = VCFBasics.VCF(line)
    qual = v.value('qual')
    if args.quality and qual < args.quality: continue
    #print line.rstrip()
    info = v.value('info')
    m0 = re.search('DP=(\d+)',info)
    dtot = int(m0.group(1))
    if args.mintotaldepth:
      if dtot < args.mintotaldepth: continue
    m = re.search('DP4=(\d+),(\d+),(\d+),(\d+)',info)
    d = [int(m.group(1)),int(m.group(2)),int(m.group(3)),int(m.group(4))]
    #if args.biallelic: 
    refhq = d[0]+d[1]
    althq = d[2]+d[3]
    if args.mindepth:
      lowest = min(refhq,althq)
      if lowest < args.mindepth: continue
    if args.mindepthstrand:
      loweststrand = min(d)
      if loweststrand < args.mindepthstrand: continue
    if args.biallelic:
      if len(v.value('ref')) != 1 or len(v.value('alt')) != 1: continue
      r = v.value('remainder').split()[1]
      m=re.match('^(\d+)[\|/](\d+):',r)
      if not m: 
        sys.stderr.write("error parsing sample field\n")
        sys.exit()
      if int(m.group(1)) == int(m.group(2)): continue
    print line.rstrip()
if __name__=="__main__":
  main()
