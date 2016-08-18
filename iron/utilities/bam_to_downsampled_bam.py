#!/usr/bin/python
import sys, argparse, gzip, os,re 
from random import shuffle
from subprocess import PIPE, Popen
from uuid import uuid4

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('bam_input',help="BAM file")
  parser.add_argument('count',type=int,help="number of read names to select")
  parser.add_argument('--aligned',action='store_true',help="only consider aligned reads")
  parser.add_argument('-o','--output',required=True,help="Output file basename, program will add .<size>.<RNG>.bam")
  args = parser.parse_args()

  id = uuid4().hex[:6]

  sys.stderr.write("Traversing bam to discover names\n")
  cmd = 'samtools view '+args.bam_input
  if args.aligned: cmd += ' -F 4'
  names = []
  p = Popen(cmd.split(),stdout=PIPE)
  prog = re.compile('^([^\t]+)')
  z = 0
  for line in p.stdout:
    z += 1
    if z%10000==0: sys.stderr.write("traversed "+str(z)+"   \r")
    m = prog.match(line)
    names.append([m.group(1),z])
  names.sort(key=lambda x: x[0])
  unique_names = list(set([x[0] for x in names]))
  shuffle(unique_names)
  sample_names = sorted(unique_names[0:args.count])
  #for name in sample_names:
  #  print name
  lines = []
  nind = 0
  sind = 0
  totalcnt = 0
  while sind < len(sample_names):
    if sample_names[sind] == names[nind][0]: 
      totalcnt += 1
      while True:
        lines.append(names[nind][1])
        if nind+1==len(names): break #already on last
        if names[nind+1][0]!=names[nind][0]: break # break on different name
        nind+=1
      sind += 1
    nind+=1
  sys.stderr.write("\n")
  sys.stderr.write("getting "+str(len(lines))+" total lines\n")
  sys.stderr.write("getting "+str(totalcnt)+" total read names\n")
  cmd = "samtools view -Sb - -o "+args.output+'.'+str(totalcnt)+'.'+id+'.bam'
  p0 = Popen(cmd.split(),stdin=PIPE,close_fds=True)
  cmd = "samtools view -H "+args.bam_input
  p1 = Popen(cmd.split(),stdout=p0.stdin)
  p1.communicate()
  cmd = 'samtools view '+args.bam_input
  if args.aligned: cmd += ' -F 4'
  p = Popen(cmd.split(),stdout=PIPE)
  z = 0
  lines.sort()
  lind = 0
  for line in p.stdout:
    z += 1
    if z%10000==0: sys.stderr.write("Final traversal of "+str(z)+" alignments    \r")
    if lind >= len(lines): continue
    if lines[lind] == z:
      p0.stdin.write(line) # write to the p0 pipe
      lind+=1
  p.communicate()
  sys.stderr.write("\n")
  p0.communicate()  
if __name__=="__main__":
  main()
