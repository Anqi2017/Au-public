#!/usr/bin/python
import argparse, sys
from GenePredBasics import GenePredEntry
import RangeBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input_gpd',help="GENEPRED input or - for STDIN")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input_gpd != '-': inf = open(args.input_gpd)
  seen = set()
  ls = RangeBasics.Loci()
  ls.verbose = True
  ls.use_direction = False
  for line in inf:
    if line[0] == '#': continue
    gpd = GenePredEntry(line)
    if gpd.value('name') in seen:
      sys.stderr.write("ERROR: need uniquely named genepred entry names\n"+name+"\n")
      sys.exit()
    seen.add(gpd.value('name'))
    r = gpd.locus_range.copy()
    r.direction = None
    r.set_payload(gpd.value('name'))
    l = RangeBasics.Locus()
    l.add_member(r)
    ls.add_locus(l)
  ls.update_loci()
  z = 0
  for locus in ls.loci:
    z += 1
    for member in locus.members:
      print str(z) + "\t" + member.get_payload()
if __name__=="__main__":
  main()
