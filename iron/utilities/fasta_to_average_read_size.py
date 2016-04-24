#!/usr/bin/python
import sys, argparse
from SequenceBasics import FastaHandleReader
#Pre: Input is
#     1. An alignment to a transcriptome in query order
#Post: Count Avg Stddev written to STDERR on the fly

def main():
  parser = argparse.ArgumentParser(description="Find mapping distance of paired end reads.  Takes an ordered (by query) alignment to a transcriptome.\nSomething that works for an input thus far is like:\nhisat --reorder -x mytranscriptome -1 my_1.fastq -2 my_2.fastq | this_script.py -")
  parser.add_argument('input_fasta',help="FASTAFILE or - for stdin")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input_fasta != '-':
    inf = open(args.input_fasta)
  fh = FastaHandleReader(inf)
  data = []
  sys.stderr.write("Reads    Mean    Stddev\n")
  while True:
    entry = fh.read_entry()
    if not entry: break
    dist = len(entry['seq'])
    data.append(dist)
    if len(data) < 2: continue
    if len(data) %1000 ==0: sys.stderr.write(str(len(data))+"    "+str(int(mean(data)))+"    "+str(int(stddev(data)))+"              \r")
  sys.stderr.write(str(len(data))+"    "+str(int(mean(data)))+"    "+str(int(stddev(data)))+"              \r")
  sys.stderr.write("\n")
def mean(data):
  n = len(data)
  if n < 1:
    sys.stderr.write("ERROR: need at least one data point\n")
  return sum(data)/float(n)

def _ss(data):
  c = mean(data)
  ss = sum([(x-c)**2 for x in data])
  return ss

def stddev(data):
  n = len(data)
  if n < 2:
    sys.stderr.write("ERROR: need at least two data points\n")
  ss = _ss(data)
  var = ss/n
  return var**0.5 

if __name__=="__main__":
  main()
