#!/usr/bin/python
import sys, argparse, re
from FileBasics import GenericFileReader

# Pre: Takes a genepred file name as an input and an output file name for where to write the output
# Post:  Outputs a genepred file to your output location that now has each transcript uniquely named
#        Also outputs a file at your outputlocation ".key_file" that contains a key
#        with unique_transcript_name followed by the original genepred line.
#  Modifies: FIle IO

def main():
  parser = argparse.ArgumentParser(description="Make a gpd with unique transcript names and a key to their original gpd entry\n")
  parser.add_argument("gpd_infile",help="FILENAME genepred file")
  parser.add_argument("gpd_outfile",help="FILENAME genepred file")
  args = parser.parse_args()
  gfr = GenericFileReader(args.gpd_infile)
  seen = {}
  while True:
    line = gfr.readline()
    if not line: break
    if re.match('^#',line): continue
    line = line.rstrip()
    f = line.split("\t")
    if f[1] not in seen:
      seen[f[1]] = []
    seen[f[1]].append(line)
  gfr.close()
  of_gpd = open(args.gpd_outfile,'w')
  of_key = open(args.gpd_outfile+".key_file",'w')
  for tx in seen:
    for i in range(0,len(seen[tx])):
      name = tx
      if len(seen[tx]) > 1:
        name = tx + '['+str(i+1)+'/'+str(len(seen[tx]))+']'
      f = seen[tx][i].split("\t")
      f[1] = name
      newline = "\t".join(f)
      of_key.write(name + "\t" + seen[tx][i] + "\n")       
      of_gpd.write(newline + "\n")       
  of_key.close()
  of_gpd.close()
main()
