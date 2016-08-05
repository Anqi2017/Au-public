#!/usr/bin/python
import sys, argparse, gzip, os.path
#from Bio.Format.Sam import BAMFile
from subprocess import Popen, PIPE

def main():
  parser = argparse.ArgumentParser(description="Get the best alignments from an indexed bam file",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use bgi (our indexing) BAM file input")
  parser.add_argument('-o','--output',help="The output file or SAM STDOUT if unset")
  args = parser.parse_args()

  if args.output:
    cmd = 'samtools view -Sb -'
    of2 = open(args.output,'w')
    p = Popen(cmd.split(),stdin=PIPE,stdout=of2)
    of = p.stdin
  else:
    of = sys.stdout

  sys.stderr.write("Checking index...\n")
  if not os.path.isfile(args.input+'.bgi'):
    sys.stderr.write("ERROR: bgi index (our format) needs to be set\n")
    sys.exit()
  # get best
  inf = gzip.open(args.input+'.bgi')
  z = 0
  bestlines = []
  for line in inf:
    z += 1
    f = line.rstrip().split("\t")
    flag = int(f[5])
    if flag & 2304 != 0:
      continue
    bestlines.append(z)

  inf.close()
  #bf = BAMFile(args.input)
  cmd = 'samtools view -H '+args.input
  pin = Popen(cmd.split(),stdout=of)
  pin.communicate()
  sys.stderr.write("Traversing bam file...\n")
  k=0
  #tot = bf.index.get_length()
  #of.write(bf.header_text)
  cmd = 'samtools view '+args.input
  pin = Popen(cmd.split(),stdout=PIPE)
  bestiter = 0
  total = len(bestlines)
  for line in pin.stdout:
    k+=1
    if k%1000==0:sys.stderr.write(str(k)+'/'+str(z)+"\r")
    if bestiter >= total: continue
    if k != bestlines[bestiter]:
      continue
    bestiter+=1
    of.write(line)
  pin.communicate()
  sys.stderr.write("\n")
  if args.output:
    p.communicate()
    of2.close()
  else:
    of.close()
  pin.stdout.close()
if __name__=="__main__":
  main()
