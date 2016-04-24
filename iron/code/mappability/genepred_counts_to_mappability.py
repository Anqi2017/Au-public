#!/usr/bin/python
import sys, argparse, re, os
from SequenceBasics import decode_name
from multiprocessing import cpu_count, Pool, Lock
from subprocess import Popen, PIPE
from SamBasics import SAM, is_header
from StatisticsBasics import average, median

null = open(os.devnull,'w')
glock = Lock()

def main():
  parser = argparse.ArgumentParser(description="For every genepred entry report its alignability",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="STDIN is -")
  parser.add_argument('-k','--fragment_size',default=100,type=int,help="Fragment size to try to align")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="number of threads")
  parser.add_argument('--type',choices=['mean','median'],default='mean',help="how to combine mappability fraction")
  parser.add_argument('--perbase',action='store_true',help='show all averages')
  parser.add_argument('--output','-o',help="output file or STDOUT if not set")
  args = parser.parse_args()
  
  if args.input == '-': args.input = sys.stdin
  else: args.input = open(args.input)

  if args.output:  args.output = open(args.output,'w')
  else: args.output = sys.stdout

  buffer = []
  prev = -1
  for line in args.input:
    if is_header(line): continue
    sam = SAM(line)
    name = decode_name(sam.value('qname')).split("\t")
    qlen = len(sam.value('seq'))
    if qlen != args.fragment_size:
      sys.stderr.write("WARNING qlen != fragment_size\n")
    cnt = 0
    if sam.value('cigar')!='*': 
      m = re.search('NH:i:(\d+)',sam.value('remainder'))
      if not m: 
        sys.stderr.write("ERROR not hisat format\n")
        sys.exit()
      cnt = int(m.group(1))
    name[2] = int(name[2])
    name[3] = int(name[3])
    name[4] = int(name[4])
    name.append(cnt)
    name.append(qlen)
    if name[2] != prev:
      if len(buffer) > 0: 
        output_buffer(buffer,args)
        buffer = []
    prev = name[2]
    buffer.append(name)
  if len(buffer) > 0:
    output_buffer(buffer,args)


def output_buffer(buffer,args):
  if len(buffer) == 0: return
  name = buffer[0][0]
  gene_name = buffer[0][1]
  line_number = buffer[0][2]
  #get the values for the positions
  positions = {}
  for e in buffer:
    origlen = e[3]
    pos = e[4]
    cnt = e[5]
    seqlen = e[6]
    if pos not in positions: positions[pos] = [0,seqlen]
    if cnt > positions[pos][0]: positions[pos] = [cnt,seqlen]
  # should have one number per position
  allcnts = {}
  #print len(positions.keys())
  for pos in sorted(positions.keys()):
    #print pos
    for i in range(pos,pos+positions[pos][1]): #for the length of this read
      if i not in allcnts: allcnts[i] = []
      allcnts[i].append(positions[pos][0]) # add the count to this position
  #now we should have all counts for all positions
  avs = []
  for pos in sorted(allcnts.keys()):
    cnts = allcnts[pos]
    nonzeros = [x for x in cnts if x != 0]
    zeros = [x for x in cnts if x == 0]
    estimatezero = cnts[:]
    for i in range(0,len(estimatezero)): 
      if estimatezero[i] == 0: estimatezero[i] = 100000000
    # classify it as unmappable if the majority of reads are unmappable
    #if len(zeros) > len(nonzeros): 
    #  avs.append(float(0))
    #  continue
    fracs = 0
    if len(estimatezero) > 0:
     if  args.type == 'median':
       fracs = median([float(1)/float(x) for x in estimatezero])
     elif args.type == 'mean':
       fracs = average([float(1)/float(x) for x in estimatezero]) 
    avs.append(float(fracs))
  #avs contains our average mappability
  zeros = [x for x in avs if x<0.01]
  multis = [x for x in avs if x >= 0.01 and x <= 0.5]
  singles = [x for x in avs if x > 0.5]
  ln = name+"\t"+gene_name+"\t"+str(origlen)+"\t"+str(len(zeros))+"\t"+str(len(multis))+"\t"+str(len(singles))
  if args.perbase:  ln += "\t"+','.join([str(round(x,4)) for x in avs])
  args.output.write(ln+"\n")
if __name__=="__main__":
  main()
