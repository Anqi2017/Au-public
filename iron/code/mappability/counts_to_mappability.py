#!/usr/bin/python
import sys, argparse, re
from math import log10, pow

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--fragment_length',type=int,required=True)
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--median',action='store_true')
  group.add_argument('--mean',action='store_true')
  group.add_argument('--geometric_mean',action='store_true')
  args = parser.parse_args()
  
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)

  k = args.fragment_length
  buffer = []
  for i in range(0,k):
    buffer.append([])
  prev_chr = ''
  for line in args.input:
    m = re.match('^([^:]+):(\d+)-(\d+)\s+(\d+)$',line.rstrip())
    if not m: 
      sys.stderr.write("ERROR unexpected read format\n")
      sys.exit()
    chr = m.group(1)
    start = int(m.group(2))-1
    end = int(m.group(3))
    cnt = int(m.group(4))
    if prev_chr != chr:
      #we have a new chromosome
      buffer = []
      for i in range(0,k): buffer.append([])
    prev_chr = chr
    prev = (start-1)%k
    buffer[prev] = []
    for i in range(start,end):
      buffer[i%k].append(cnt) 
    curr = buffer[start%k]
    nonzero = [x for x in curr if x != 0]
    zero = [x for x in curr if x == 0]
    value = 0
    # require at least half the reads to be nonzero
    if len(nonzero) >= len(zero): value = custom_avg(nonzero,args)
    print chr +"\t"+str(start)+"\t"+str(start+1)+"\t"+str(value)

def custom_avg(nonzero,args):
  fracs = [float(1)/float(x) for x in nonzero]
  if args.median:  return median(fracs)
  elif args.mean: return float(sum(fracs))/float(len(fracs))
  elif args.geometric_mean: return pow(10,float(sum([log10(float(x)) for x in fracs]))/float(len(fracs)))
  else:
    sys.stderr.write("ERROR unknown args\n")
    sys.exit()
    
def median(arr):
  if len(arr) == 1: return arr[0]
  quot = len(arr)/2
  rem = len(arr)%2
  if rem != 0:
    return sorted(arr)[quot]
  return float(sum(sorted(arr)[quot-1:quot+1]))/float(2)
if __name__=="__main__":
  main()
