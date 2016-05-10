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

  bf = BAMFile(args.input,skip_index=True)
  sys.stderr.write("Checking index...\n")
  if not os.path.isfile(args.input+'.bgi'):
    sys.stderr.write("ERROR: bgi index (our format) needs to be set\n")
    sys.exit()
  inf = gzip.open(args.input+'.bgi')
  k = 0
  best = {}
  for line in inf:
    k += 1
    f = line.rstrip().split("\t")
    if f[0] not in best:
      best[f[0]] = [k,int(f[4])]
    if f[4] > best[f[0]][1]:
      best[f[0]] = [k,int(f[4])]
  tot = k
  inf.close()
  sys.stderr.write("Traversing bam file...\n")
  bestlines = {}
  for name in best:  bestlines[best[name][0]] = name
  args.output.write(bf.header_text.rstrip()+"\n")
  k = 0
  for e in bf:
    k+=1
    if k%1000==0:sys.stderr.write(str(k)+'/'+str(tot)+"\r")
    if k in bestlines:
      if e.value('qname') != bestlines[k]:
        sys.stderr.write("\nERROR not a 1:1 expected correspondence between index lines and entries\n")
      args.output.write(e.get_line().rstrip()+"\n")
  sys.stderr.write("\n")
if __name__=="__main__":
  main()
