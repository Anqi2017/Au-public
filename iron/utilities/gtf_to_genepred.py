#!/usr/bin/python
import GTFBasics
import sys

# Pre:  A GTF filename, for a file with 'exon' features, and 'gene_id' and 'transcript_id' attributes.
# Post: Prints to stdout a genepred with transcripts

def main():
  if len(sys.argv) < 2: 
    sys.stderr.write("gtf_to_genepred.py <gtf filename>\n")
    return
  gtf = GTFBasics.GTFFile(sys.argv[1])
  gtf.write_genepred(sys.stdout)

main()
