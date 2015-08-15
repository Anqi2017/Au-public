#!/usr/bin/python
import argparse, sys
import SamBasics

def main():
  parser = argparse.ArgumentParser(description="Convert a sam file into a psl file")
  parser.add_argument('--genome',help="FASTA input file of reference genome")
  parser.add_argument('infile',help="FILENAME input file or '-' for STDIN")
  args = parser.parse_args()
  inf = sys.stdin
  if args.infile != '-': 
    inf = open(args.infile)
  spcf = SamBasics.SAMtoPSLconversionFactory()
  if args.genome: spcf.set_genome(args.genome)
  for line in inf:
    line = line.rstrip()
    if SamBasics.is_header(line): 
      spcf.read_header_line(line)
      continue
    psl = spcf.convert_line(line)
    if psl:
      print psl
  inf.close()


if __name__=="__main__":
  main()
