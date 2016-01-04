#!/usr/bin/python
import sys, argparse
import PSLBasics
from SequenceBasics import read_fasta_into_hash, FastaHandleReader, rc
# Fixes the matches/mismatches/Ncount of a PSL file
# This file will take as an input, a PSL file,
# a query file, and a reference file.
# Pre: Every query sequence is required to be present
#      in the query file and every reference in the reference file
# IMPORTANT!!!! input PSL is ordered by query, and query FASTA is ordered by query

def main():
  parser = argparse.ArgumentParser(description="Correct the matches/mismatches and Ncount of a PSL file")
  parser.add_argument('input',help="PSLFILE or - for STIDN")
  parser.add_argument('reference',help="FASTAFILE reference genome")
  parser.add_argument('query',help="FASTAFILE query sequences")
  parser.add_argument('--minimum_intron_size',type=int,default=68,help="INT")
  #parser.add_argument('--ordered_query',action='store_true',help="The query psl and fasta are both ordered by query name for optimal performance")
  args = parser.parse_args()
  # Read in the reference genome
  sys.stderr.write("Reading in reference genome\n")
  g = read_fasta_into_hash(args.reference)
  sys.stderr.write("Finished reading "+str(len(g.keys()))+" reference sequences\n")
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  fhr = FastaHandleReader(open(args.query))
  last_fasta = fhr.read_entry()
  if not last_fasta:
    sys.stderr.write("ERROR: No query sequences\n")
    sys.exit()
  for line in inf:
    p = PSLBasics.PSL(line)
    if not p.validate():
      sys.stderr.write("WARNING skipping invalid PSL entry. This script fixes improper mismatch and match counts .. and gap counts... it doesn't perform miracles. Problem line:\n"+line.rstrip()+"\n")
    n = p.value('qName')
    if not last_fasta:
      sys.stderr.write("ERROR: Ran out of query sequences too soon.  Are they sorted properly\n")
      sys.exit()
    while last_fasta['name'] != n:
      last_fasta = fhr.read_entry()
    p.set_query(last_fasta['seq'])
    p.set_reference_dictionary(g)
    print p.get_line()
    p.pretty_print(50)
  fhr.close()
if __name__=="__main__":
  main()
