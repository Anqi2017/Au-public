#!/usr/bin/python
import argparse, sys
from SamBasics import is_header
from subprocess import Popen, PIPE

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="input sorted bam or - for STDIN. expects header")
  parser.add_argument('--positional_duplicates',type=int,help="maximum number of positional duplicatse to allow through from a sorted sam")
  args = parser.parse_args()
  in_header = True
  bam = False
  if args.input == '-': args.input = sys.stdin
  else: 
    bam = True
    cmd = "samtools view -h "+args.input
    p = Popen(cmd.split(),stdout=PIPE)
    args.input = p.stdout
  line = args.input.readline()
  if not line: return #done
  buffer_name = ''
  buffer_count = 0
  while True:
    if in_header: 
        line = args.input.readline()
        if not line: break
        if is_header(line):
          print line.rstrip()
          continue
        else: in_header = False
    else:
      line = args.input.readline()
      if not line: break
    #have a line
    f = line.split("\t")
    if args.positional_duplicates:
      pos = ':'.join([f[2],f[3],f[5]])
      if pos != buffer_name:
        buffer_count = 0
        buffer_name = pos
      buffer_count += 1
      if buffer_count <= args.positional_duplicates:
        print line.rstrip()
    else:
      print line.rstrip()
  if bam: p.communicate()

if __name__=="__main__":
  main()
