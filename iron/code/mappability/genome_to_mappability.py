#!/usr/bin/python
import sys, argparse, os
from SequenceBasics import FastaHandleReader
from subprocess import Popen, PIPE
from multiprocessing import cpu_count

def main():
  parser = argparse.ArgumentParser(description="Convert a genome to its mappability",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('reference_genome',help="Use - for STDIN")
  parser.add_argument('-k','--fragment_length',type=int,default=36,help="length of fragment to check mappability")
  parser.add_argument('-x','--genome_index',required=True)
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Thread count")
  parser.add_argument('-o','--output',help="set for output file otherwise will be STDOUT")
  parser.add_argument('--type',choices=['mean','median','geometric_mean'],default='mean',help="How to combine window results")
  args = parser.parse_args()

  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout

  
  
  udir = os.path.dirname(os.path.realpath(__file__))
  cmd4 = 'bed_tools.py - --merge --break_merge_on_feature'
  p4 = Popen(cmd4.split(),stdin=PIPE,stdout=args.output)
  cmd3 = udir+'/counts_to_mappability.py - --fragment_length '+str(args.fragment_length)
  cmd3 += ' --'+args.type
  p3 = Popen(cmd3.split(),stdin=PIPE,stdout=p4.stdin)
  cmd2 = 'hisat_to_mapping_count.py -'
  p2 = Popen(cmd2.split(),stdin=PIPE,stdout=p3.stdin)
  cmd1 = 'hisat -x '+args.genome_index+' -U - -f --reorder -p '+str(args.threads)
  p1 = Popen(cmd1.split(),stdin=PIPE,stdout=p2.stdin)
  inf = open(args.reference_genome)
  fhr = FastaHandleReader(inf)
  while True:
    e = fhr.read_entry()
    if not e: break
    for i in range(0,len(e['seq'])-args.fragment_length):
      p1.stdin.write('>'+e['name']+':'+str(i+1)+'-'+str(i+args.fragment_length)+"\n")
      p1.stdin.write(e['seq'][i:i+args.fragment_length].upper()+"\n")
  p1.communicate()
  p2.communicate()
  p3.communicate()
  p4.communicate()
  args.output.close()

if __name__=="__main__":
  main()
