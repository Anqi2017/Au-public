#!/usr/bin/python
import sys, argparse, re, os
from subprocess import Popen, PIPE
from multiprocessing import Process

def main():
  parser = argparse.ArgumentParser(description="Split aligned reads into chromosomes",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="sorted samtools indexed BAM file")
  parser.add_argument('-o','--output',help="output directory")
  args = parser.parse_args()
  m = re.search('([^\/]+)\.bam$',args.input)
  if not m: 
    sys.stderr.write("bam file must end in .bam")
    sys.exit()
  basename = m.group(1)
  args.output = args.output.rstrip('/')
  if not os.path.exists(args.output):
    os.makedirs(args.output)
  cmd = "samtools view -H "+args.input
  p = Popen(cmd.split(),stdout=PIPE)
  chrs = []
  for line in p.stdout:
    m = re.match('@SQ\tSN:(\S+)',line)
    if not m: continue
    chrs.append(m.group(1))
  p.communicate()
  ps = []
  for chr in chrs:
    ps.append(Process(target=do_chr,args=(chr,args,basename,)))
  ps.append(Process(target=do_un,args=(args,basename,)))
  for p in ps: p.start()
  for p in ps: p.join()

def do_un(args,basename):
    of = open(args.output+'/'+basename+'.unaligned.bam','w')
    cmd1 = "samtools view -f 4 -h "+args.input
    p1 = Popen(cmd1.split(),stdout=PIPE)
    cmd2 = "samtools view -Sb -"
    p2 = Popen(cmd2.split(),stdin=p1.stdout,stdout=of)
    p2.communicate()
    p1.communicate()
    of.close()

def do_chr(chr,args,basename):
    of = open(args.output+'/'+basename+'.'+chr+'.bam','w')
    cmd1 = "samtools view -F 4 -h "+args.input+" "+chr
    p1 = Popen(cmd1.split(),stdout=PIPE)
    cmd2 = "samtools view -Sb -"
    p2 = Popen(cmd2.split(),stdin=p1.stdout,stdout=of)
    p2.communicate()
    p1.communicate()
    of.close()

if __name__=="__main__":
  main()
