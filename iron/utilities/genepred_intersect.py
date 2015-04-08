#!/usr/bin/python
import GenePredBasics, RangeBasics
import argparse, sys

# This should be a powerful commandline utility understanding the overlaps of genepred files

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-a',required=True,help='FILENAME genepred file A')
  parser.add_argument('-b',required=True,help='FILENAME genepred file B')
  parser.add_argument('--minexoncount',nargs='?',help='INT the minimum number of exons required.')
  parser.add_argument('--leftouterjoin',nargs='?',help='OUTPUT the entry A regardless of whether a matching entry in B is found')
  parser.add_argument('--minoverlap_internal',nargs='?',help='FLOAT the fraction (0-1) of the required reciprocal overlap of an internal exon to call an exon a match.')
  parser.add_argument('--minoverlap_first',nargs='?',help='FLOAT the fraction (0-1) of the required reciprocal overlap of the first exon to call an exon a match.')
  parser.add_argument('--minoverlap_last',nargs='?',help='FLOAT the fraction (0-1) of the required reciprocal overlap of the last exon to call an exon a match.')
  parser.add_argument('--minoverlap',nargs='?',help='FLOAT the fraction (0-1) of the required reciprocal overlap of any exon to call an exon a match.')
  parser.add_argument('--output_a_not_in_b',action='store_true',help='Output entries that occur in A but not B')
  parser.add_argument('--best_b_only',action='store_true',help='Output only one entry of B for each A and try to pick the best based on reciprocal overlap')
  parser.add_argument('--allow_a_subset_of_b_fragments',action='store_true',help='If A is just a subset of B, then call it as a match.')
  parser.add_argument('--allow_any_fragments',action='store_true',help='If set, allow any partial match, not just the best')
  args = parser.parse_args()
  
  # go through contingencies of overlap requirements and set them
  overlap = [0,0,0]
  if args.minoverlap:
    overlap = [float(args.minoverlap), float(args.minoverlap), float(args.minoverlap)]
  if args.minoverlap_first:
    overlap[0] = float(args.minoverlap_last)
  if args.minoverlap_last:
    overlap[2] = float(args.minoverlap_last)
  if args.minoverlap_internal:
    overlap[1] = float(args.minoverlap_internal)

  # read the genepred files
  gpdA = GenePredBasics.GenePredFile(args.a)
  gpdB = GenePredBasics.GenePredFile(args.b)

  for eA in gpdA.entries:
    a_unique = True
    best_exon_count = 0
    best_overlap = 0
    best_line = ''
    for eB in gpdB.entries:
      double_line = GenePredBasics.entry_to_line(eA.entry) + "\t" + GenePredBasics.entry_to_line(eB.entry) + "\n"
      gpd_comparison = GenePredBasics.GenePredComparison()
      gpd_comparison.set_overlap_requirement(overlap)
      # normal is to do full length matches
      if not (args.allow_a_subset_of_b_fragments or args.allow_any_fragments):
        gpd_comparison.set_require_all_exons_overlap(True)
        gpd_comparison.compare(eA,eB)
        if gpd_comparison.output['full_match']:
          a_unique = False
          if args.output_a_not_in_b:
            break # we can bust out of the inner loop if we are only printing stuff unique to a 
          if not args.best_b_only: # if we aren't waiting for the best, print it
            sys.stdout.write(double_line)
          else:
            # only do the best
            if gpd_comparison.output['consecutive_exons'] > best_exon_count \
            or (gpd_comparison.output['consecutive_exons'] == best_exon_count \
            and gpd_comparison.output['overlap_length'] > best_overlap):
              best_exon_count = gpd_comparison.output['consecutive_exons']
              best_overlap = gpd_comparison.output['overlap_length']
              best_line = double_line

      # Allow partial matches
      else:          
        gpd_comparison.compare(eA,eB)
        if gpd_comparison.output['partial_match']:
          # if we require a to be subset of b
          if args.allow_a_subset_of_b_fragments \
          and not (eA.get_exon_count() < eB.get_exon_count() \
          and eA.get_exon_count() == gpd_comparison.output['consecutive_exons']):
            break
          a_unique = False
          if args.output_a_not_in_b:
            break
            # only do the best
          if not args.best_b_only:
            sys.stdout.write(double_line)
          else:
            if gpd_comparison.output['consecutive_exons'] > best_exon_count \
            or (gpd_comparison.output['consecutive_exons'] == best_exon_count \
            and gpd_comparison.output['overlap_length'] > best_overlap):
              best_exon_count = gpd_comparison.output['consecutive_exons']
              best_overlap = gpd_comparison.output['overlap_length']
              best_line = double_line
    if best_exon_count > 0 and args.best_b_only:
      sys.stdout.write(best_line)
    if a_unique and (args.output_a_not_in_b or args.leftouterjoin):
      sys.stdout.write(GenePredBasics.entry_to_line(eA.entry)+"\n")

def compare(eA,eB,overlap):
  range_list_A = eA.range_set.get_range_list()
  first_exon_A = eA.get_first_exon_genomic_range()
  last_exon_A = eA.get_last_exon_genomic_range()

  range_list_B = eB.range_set.get_range_list()
  first_exon_B = eB.get_first_exon_genomic_range()
  last_exon_B = eB.get_last_exon_genomic_range()
  for r in range_list_A:
    over = eB.range_set.get_overlapped(r)

    reqover = 0
    if r.equals(first_exon_B): reqover = overlap[0]
    elif r.equals(last_exon_B): reqover = overlap[2]
    else: reqover = overlap[1]

    if len(over.get_range_list()) != 1: return [0, None] # if we don't have a 1:1

    # Check our reciprocal overlap here
    osize_bp = over.get_range_list()[0].overlap_size(r)
    size1 = over.get_range_list()[0].length()
    size2 = r.length()
    if size1 == 0 or size2 ==0 or osize_bp ==0: 
      sys.stderr.write("warning strange size 0 for a passing match\n")
      return [0, None]
    longest = sorted([size1,size2])[1]
    obs_overlap = float(osize_bp)/float(longest)
    if obs_overlap < reqover:
      return [0, None]

  for r in range_list_B:
    over = eA.range_set.get_overlapped(r)
    if len(over.get_range_list()) != 1: return [0, None]

    # Check our reciprocal overlap here
    osize_bp = over.get_range_list()[0].overlap_size(r)
    size1 = over.get_range_list()[0].length()
    size2 = r.length()
    if size1 == 0 or size2 ==0 or osize_bp ==0: 
      sys.stderr.write("warning strange size 0 for a passing match\n")
      return [0, None]
    longest = sorted([size1,size2])[1]
    obs_overlap = float(osize_bp)/float(longest)
    if obs_overlap < reqover:
      return [0, None]

  return True
  
main()
