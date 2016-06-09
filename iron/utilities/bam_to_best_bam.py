#!/usr/bin/python
import sys, argparse, gzip, os.path
from Bio.Format.Sam import BAMFile

def main():
  parser = argparse.ArgumentParser(description="Get the best alignments from an indexed bam file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use bgi (our indexing) BAM file input")
  parser.add_argument('-o','--output',help="The output file or STDOUT if unset")
  args = parser.parse_args()

  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout  

  sys.stderr.write("Checking index...\n")
  if not os.path.isfile(args.input+'.bgi'):
    sys.stderr.write("ERROR: bgi index (our format) needs to be set\n")
    sys.exit()
  bf = BAMFile(args.input)
  bf.read_index()
  sys.stderr.write("Traversing bam file...\n")
  k=0
  tot = bf.index.get_length()
  args.output.write(bf.header_text)
  for e in bf:
    k+=1
    if k%1000==0:sys.stderr.write(str(k)+'/'+str(tot)+"\r")
    if not e.indexed_as_primary_alignment(): continue
    args.output.write(e.get_line().rstrip()+"\n")
  sys.stderr.write("\n")
if __name__=="__main__":
  main()
