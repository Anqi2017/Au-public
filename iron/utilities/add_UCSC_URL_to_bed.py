#!/usr/bin/python
import UCSCBasics

import argparse, sys

def main():
  parser = argparse.ArgumentParser(description='Add a UCSC URL to a bed file')
  parser.add_argument('--session',help="STRING hgsid=XXXXXXXXXXX where X's are your string")
  parser.add_argument('--first',action='store_true',help="Add it to the beginning of the file instead of the end")
  parser.add_argument('--excel',action='store_true',help="Make as an excel formated hyperlink")
  parser.add_argument('--db',help="DATABASE like hg19",required=True)
  parser.add_argument('bed_file',help="FILENAME if - then STDIN")
  args = parser.parse_args()

  urlfactory = UCSCBasics.URLfactory(args.db)
  if args.session: urlfactory.set_session(args.session)
  bed_handle = sys.stdin
  if args.bed_file != '-': bed_handle = open(args.bed_file)
  for line in bed_handle:
    f = line.rstrip().split("\t")
    chr = f[0]
    start = int(f[1])
    stop = int(f[2])
    url = urlfactory.url_from_bed_coordinates(chr,start,stop)
    if args.excel:
      url = urlfactory.excel_from_bed_coordinates(chr,start,stop)
    if args.first:
      print url + "\t" + line.rstrip()
    else:
      print line.rstrip() + "\t" + url

if __name__=="__main__":
  main()
