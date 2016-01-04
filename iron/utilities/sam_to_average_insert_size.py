#!/usr/bin/python
import sys, argparse
import SamBasics, PSLBasics

#Pre: Input is
#     1. An alignment to a transcriptome in query order
#Post: Count Avg Stddev written to STDERR on the fly

def main():
  parser = argparse.ArgumentParser(description="Find mapping distance of paired end reads.  Takes an ordered (by query) alignment to a transcriptome.\nSomething that works for an input thus far is like:\nhisat --reorder -x mytranscriptome -1 my_1.fastq -2 my_2.fastq | this_script.py -")
  parser.add_argument('input_sam',help="SAMFILE ordered alignment a transcriptome or - for stdin")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input_sam != '-':
    inf = open(args.input_sam)
  msr = SamBasics.MultiEntrySamReader(inf)
  spcf = SamBasics.SAMtoPSLconversionFactory()
  data = []
  sys.stderr.write("Pairs    Mean    Stddev\n")
  while True:
    entries = msr.read_entries()
    if not entries: break
    if len(entries)!=2: continue
    [e1,e2] =entries
    if e1.check_flag(4) or e2.check_flag(4): continue
    if not e1.check_flag(2) and e2.check_flag(2): continue
    if not ((e1.check_flag(64) and e2.check_flag(128)) or (e1.check_flag(128) and e2.check_flag(64))): continue
    p1 = spcf.convert_line(e1.get_line())
    p2 = spcf.convert_line(e2.get_line())
    if not p1 or not p2: continue
    p1 = PSLBasics.PSL(p1)
    p2 = PSLBasics.PSL(p2)
    dist = max(p2.value('tEnd')-p1.value('tStart'),p1.value('tEnd')-p2.value('tStart'))
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
