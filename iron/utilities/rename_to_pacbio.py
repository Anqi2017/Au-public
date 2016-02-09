#!/usr/bin/python
import sys,argparse
from SequenceBasics import FastaHandleReader, FastqHandleReader

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="Use - for STDIN")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--fasta',action='store_true')
  group.add_argument('--fastq',action='store_true')
  args = parser.parse_args()

  if args.input=='-': args.input = sys.stdin
  else: args.input= open(args.input)

  if args.fasta:
    args.input = FastaHandleReader(args.input)
  elif args.fastq:
    args.input = FastqHandleReader(args.input)
  z = 0
  while True:
    e = args.input.read_entry()
    if not e: break
    z+=1
    name = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
    if args.fastq:
      print '@'+name
      print e['seq']
      print '+'
      print e['qual']
    elif args.fasta:
      print '>'+name
      print e['seq']
if __name__=="__main__":
  main()
