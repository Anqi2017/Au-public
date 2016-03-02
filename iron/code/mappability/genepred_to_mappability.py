#!/usr/bin/python
import sys, argparse, re, os
from SequenceBasics import read_fasta_into_hash, encode_name
from GenePredBasics import GenePredEntry as GPD
from subprocess import Popen, PIPE
from StatisticsBasics import average, median
from multiprocessing import cpu_count
import gzip

null = open(os.devnull,'w')

def main():
  parser = argparse.ArgumentParser(description="For every genepred entry report its alignability",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Genepred can be gzipped or - for STDIN")
  parser.add_argument('-r','--reference',required=True,help="Reference fasta")
  parser.add_argument('-k','--fragment_size',default=100,type=int,help="Fragment size to try to align")
  parser.add_argument('-x','--hisat_index',required=True,help="HISAT index base name")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="number of threads")
  parser.add_argument('--type',choices=['mean','median'],default='mean',help="How to bring together overlapping reads")
  parser.add_argument('--perbase',action='store_true')
  args = parser.parse_args()
  
  if args.input=='-': args.input=sys.stdin
  elif re.search('\.gz$',args.input):
    args.input = gzip.open(args.input)
  else: args.input = open(args.input)

  udir = os.path.dirname(os.path.realpath(__file__))
  cmd2 = udir+'/genepred_counts_to_mappability.py -'
  cmd2 += ' --threads '+str(args.threads)
  cmd2 += ' -k '+str(args.fragment_size)
  if args.perbase: cmd2 += ' --perbase'
  p2 = Popen(cmd2.split(),stdin=PIPE)
  ref = read_fasta_into_hash(args.reference)
  cmd1 = 'hisat -x '+args.hisat_index+' -U - -f --reorder -p '+str(args.threads)
  p1 = Popen(cmd1.split(),stdin=PIPE,stdout=p2.stdin,stderr=null)
  line_number = 0
  for line in args.input:
    line_number +=1
    gpd = GPD(line.rstrip())
    seq = gpd.get_sequence(ref)
    for i in range(0,len(seq)-args.fragment_size+1):
      info = gpd.value('name')+"\t"+gpd.value('gene_name')+"\t"+str(line_number)+"\t"+str(len(seq))+"\t"+str(i)
      einfo = encode_name(info)
      p1.stdin.write('>'+einfo+"\n")
      p1.stdin.write(seq[i:i+args.fragment_size]+"\n")
  p1.communicate()
  p2.communicate()

if __name__=="__main__":
  main()
