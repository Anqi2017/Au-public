#!/usr/bin/python
import sys, argparse
import SamBasics

def main():
  parser = argparse.ArgumentParser(description="Make sam file compatible with tools counting on a splicemap format sam file.")
  parser.add_argument('in_sam',help="FILENAME of sam file, or '-' for STDIN")
  args = parser.parse_args()
  inf = sys.stdin
  if args.in_sam != '-':
    inf = open(args.in_sam)
  for line in inf:
    line = line.rstrip()
    if SamBasics.is_header(line):
      print line
      continue
    f = line.rstrip().split("\t")
    e = SamBasics.sam_line_to_dictionary(line)
    if SamBasics.check_flag(e['flag'],4):
      continue # skip the unmapped reads
    if SamBasics.check_flag(e['flag'],16):
      f[1] = "16"
    else:
      f[1] = "0"
    f[4] = "0"
    f[6] = "*"
    f[7] = "0"
    f[8] = "0"
    print "\t".join(f)
main()
