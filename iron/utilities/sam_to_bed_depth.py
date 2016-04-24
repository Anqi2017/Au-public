#!/usr/bin/python
import sys, argparse
from SamBasics import SamLocusStream, SAMtoPSLconversionFactory
from GenePredBasics import GenePredEntry
from PSLBasics import PSL
from subprocess import Popen, PIPE
from RangeBasics import Bed

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--min_intron',type=int,default=68)
  parser.add_argument('--min_depth',type=int,default=1)
  parser.add_argument('-v','--verbose',action='store_true')
  parser.add_argument('input',help="BAM file or stream sam input with header with -")
  args = parser.parse_args()
  if args.input == '-':  inf = sys.stdin
  else:
    cmd = "samtools view -h "+args.input
    p = Popen(cmd.split(),stdout=PIPE,bufsize=1)
    inf = p.stdout
  s = SamLocusStream(inf)
  #s.junctions_only = False
  while True:
    locus = s.read_locus()
    if not locus: break
    #print locus[0].get_range_string()+"\t"+str(len(locus[1]))
    process_locus(locus[1],args)
  if args.input != '-': p.communicate()

def process_locus(locus, args):
  depth = {}
  s2psl = SAMtoPSLconversionFactory()
  unique = {}
  chr = locus[0].value('rname')
  for sam in locus:
    p = PSL(s2psl.convert_line(sam.get_line()))
    g = GenePredEntry(p.get_genepred_line())
    g = g.get_smoothed(args.min_intron)
    for i in range(0,g.get_exon_count()):
      rng = str(g.value('exonStarts')[i])+"\t"+str(g.value('exonEnds')[i])
      if rng not in unique: unique[rng] = 0
      unique[rng]+=1
  for bstr in unique:
    [start,end] = bstr.split("\t")
    for i in range(int(start),int(end)):
      if i not in depth:  depth[i] = 0
      depth[i] += unique[bstr] # add the number of these to the depth
  #now we can print the depth
  prevdepth = 0
  prevstart = None
  lasti = None
  for i in sorted(depth.keys()):
    if depth[i] < args.min_depth: continue
    if depth[i] != prevdepth: #output what we have so far if we have something
      if prevstart: 
        output_depth(chr+"\t"+str(prevstart)+"\t"+str(lasti+1)+"\t"+str(prevdepth),args)
      prevstart = i
    prevdepth = depth[i]
    lasti = i
  if prevstart:
    output_depth(chr+"\t"+str(prevstart)+"\t"+str(lasti+1)+"\t"+str(prevdepth),args)

def output_depth(ostr,args):
  if args.verbose: sys.stderr.write(" ".join(ostr.split("\t")[0:3])+"                        \r")
  print ostr

if __name__=="__main__":
  main()
