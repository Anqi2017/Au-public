#!/usr/bin/python
import argparse, sys, re
import SamBasics

def main():
  parser = argparse.ArgumentParser(description="Convert a sam file into a psl file")
  parser.add_argument('--genome',help="FASTA input file of reference genome")
  parser.add_argument('--get_secondary_alignments',action='store_true',help="Report SA:Z secondary alignments as well")
  parser.add_argument('--get_alternative_alignments',action='store_true',help="Report XA:Z alternative alignments as well")
  parser.add_argument('--get_all_alignments',action='store_true',help="Report SA:Z and XA:Z alternative alignments as well")
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
    # We have a line to convert
    psl = spcf.convert_line(line)
    if psl:
      print psl
    # Lets look for secondary alignments to convert
    if args.get_secondary_alignments or args.get_all_alignments:
      secondary_alignments = SamBasics.get_secondary_alignments(line.rstrip())
      for samline in secondary_alignments:
        psl = spcf.convert_line(samline)
        if psl:
          #print "\nsecondary"
          #print samline
          print psl
    if args.get_alternative_alignments or args.get_all_alignments:
      alternative_alignments = SamBasics.get_alternative_alignments(line.rstrip())
      for samline in alternative_alignments:
        psl = spcf.convert_line(samline)
        if psl:
          #print "\nsecondary"
          #print samline
          print psl
  inf.close()


if __name__=="__main__":
  main()
