#!/usr/bin/python
import sys, argparse, re
from Bio.Format.Sam import SamStream
from subprocess import Popen, PIPE
from multiprocessing import cpu_count

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use bam file")
  parser.add_argument('output',help="Use bam file")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Thread count")
  args = parser.parse_args()
  m = re.match('^(\S+)\.bam$',args.output)
  if not m:
    sys.stderr.write("use bam output")
    sys.exit()
  cmd1 = 'samtools view -H '+args.input
  p1 = Popen(cmd1.split(),stdout=PIPE)
  bs = SamStream(p1.stdout)
  rlens = bs.get_header().get_sequence_lengths()
  htext = bs.header_text  
  p1.communicate()
  hlines = htext.rstrip().split("\n")
  done_lens = False
  cmd = 'samtools sort -@ '+str(args.threads)+'  - '+m.group(1)
  sys.stderr.write(cmd+"\n")
  p = Popen(cmd.split(),stdin=PIPE)
  for ln in hlines:
    if re.match('@SQ\tSN:',ln):
      if not done_lens:
        done_lens = True
        for chr in sorted(rlens.keys()):
          p.stdin.write("@SQ\tSN:"+chr+"\tLN:"+str(rlens[chr])+"\n")
    else:
      p.stdin.write(ln.rstrip("\n")+"\n")
  cmd1 = 'samtools view '+args.input
  p1 = Popen(cmd1.split(),stdout=p.stdin)
  p1.communicate()
  p.communicate()

if __name__=="__main__":
  main()
