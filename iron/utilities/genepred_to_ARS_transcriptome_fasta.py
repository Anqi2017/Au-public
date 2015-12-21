#!/usr/bin/python
import sys, argparse
from ArtificalReferenceSequenceBasics import ARS, decode_ars_name
from SequenceBasics import read_fasta_into_hash
import GenePredBasics
from RangeBasics import Bed

def main():
  parser = argparse.ArgumentParser(description='Create artifical reference sequences from a genepred')
  parser.add_argument('gpd_file')
  parser.add_argument('reference_fasta')
  parser.add_argument('-o','--output',help="output file to write to or STDOUT if not set")
  args = parser.parse_args()
  of  = sys.stdout
  if args.output: of = open(args.output,'w')
  f = read_fasta_into_hash(args.reference_fasta)
  with open(args.gpd_file) as inf:
    for line in inf:
      gpd = GenePredBasics.GenePredEntry()
      gpd.line_to_entry(line.rstrip())
      ars = ARS()
      beds = []
      for i in range(0,gpd.value('exonCount')):
        b = Bed(gpd.value('chrom'),gpd.value('exonStarts')[i],gpd.value('exonEnds')[i],gpd.value('strand'))
        beds.append(b)
      ars.set_bounds(beds)
      ars.set_name(gpd.value('name'))
      ars.set_sequence_from_original_reference_hash(f)
      of.write(ars.get_fasta())

if __name__=="__main__":
  main()
