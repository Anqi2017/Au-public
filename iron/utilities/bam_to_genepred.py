#!/usr/bin/python
import sys, argparse, gzip
from Bio.Format.Sam import BAMFile, SamStream

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="BAM file or Use - for STDIN for SAM")
  parser.add_argument('--minimum_intron',type=int,default=68,help="smallest intron")
  parser.add_argument('-o','--output',help="Output file, gzip is okay")
  args = parser.parse_args()

  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz':
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')  

  if args.input =='-':
    sh = SamStream(sys.stdin)
  else:
    sh = BAMFile(args.input)
  for e in sh:
    if not e.is_aligned(): continue
    gpd_line = e.get_target_transcript(min_intron=args.minimum_intron).get_gpd_line()
    of.write(gpd_line+"\n")
  sh.close()
  of.close()
if __name__=="__main__":
  main()
