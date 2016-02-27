#!/usr/bin/python
import sys, argparse

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  group1 = parser.add_mutually_exclusive_group()
  group1.add_argument('--pad',type=int,help="add bases to bed")
  group1.add_argument('--merge',action='store_true',help="merge sorted bed entries")
  args = parser.parse_args()
  
  
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)

  # Do pads
  if args.pad:
    do_pad(args)
    return

  # Do merge
  if args.merge:
    do_merge(args)
    return

def do_merge(args):
  prev = args.input.readline()
  if not prev: return
  f = prev.rstrip().split("\t")
  buffer_chr = f[0]
  buffer_start = int(f[1])
  buffer_end = int(f[2])
  while True:
    curr = args.input.readline()
    if not curr:
      print buffer_chr+"\t"+str(buffer_start)+"\t"+str(buffer_end)
      return
    f2 = curr.rstrip().split("\t")
    #sys.stderr.write(str(f2)+"\n")
    curr_chr = f2[0]
    curr_start = int(f2[1])
    curr_end = int(f2[2])
    if curr_start > buffer_end or curr_chr !=  buffer_chr:
      #output buffer
      print buffer_chr+"\t"+str(buffer_start)+"\t"+str(buffer_end)
      buffer_chr = curr_chr
      buffer_start = curr_start
      buffer_end = curr_end
    if curr_end > buffer_end: buffer_end = curr_end
    prev = curr



def do_pad(args):
  for line in args.input:
    f = line.rstrip().split("\t")
    ln = f[0]+"\t"+str(max(0,int(f[1])-args.pad))+"\t"+str(int(f[2])+args.pad)
    if len(f) > 3:
      ln += "\t"+"\t".join(f[3:])
    print ln
  return

if __name__=="__main__":
  main()
