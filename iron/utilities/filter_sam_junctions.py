#!/usr/bin/python
import sys, argparse, re
from subprocess import Popen, PIPE
from Bio.Format.Sam import SamStream
from Bio.Stream import LocusStream
from Bio.Range import Bed
from Bio.Structure import Junction

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN or read bam")
  parser.add_argument('--minimum_intron_size',type=int,default=68,help="Require intron to be this size or larger")
  parser.add_argument('--minimum_overhang',type=int,default=10,help="At least this many bases on each side of an intron")
  parser.add_argument('--minimum_support',type=int,default=1,help="Minimum number of reads that should support a junction to report any of the reads")
  args = parser.parse_args()

  inf = sys.stdin
  if args.input != '-':
    cmd = 'samtools view -F 4 -h '+args.input
    p = Popen(cmd.split(),stdout=PIPE,bufsize=1)
    inf = p.stdout
  cmd2 = 'awk '+"'"+'{if(NF<10) print $0; else if($6~/N/) print $0;}'+"'"
  p2 = Popen(cmd2,stdout=PIPE,stdin=inf,bufsize=1,shell=True)
  stream = SamStream(p2.stdout,minimum_intron_size=args.minimum_intron_size,minimum_overhang=args.minimum_overhang)
  lstream = LocusStream(stream)
  for h in stream.header:
    print h.rstrip()
  for r in lstream:
    # now we have all the possible junctions from the range
    [juncs,sams] = get_junctions(r.get_payload(),args)
    evidence = {}
    lines = {}
    for x in set([x[0].get_range_string() for x in juncs]):
      evidence[x] = 0
      lines[x] = set()
    for i in range(0,len(juncs)):
      jstr = juncs[i][0].get_range_string()
      evidence[jstr]+=1
      lines[jstr].add(juncs[i][1])
    accepted = set()
    for jstr in evidence:
      if evidence[jstr] >= args.minimum_support:
        #print jstr
        for i in lines[jstr]: accepted.add(i)
    for i in sorted(list(accepted)):
      print sams[i].get_line().rstrip()
  p2.communicate()
  if args.input != '-':
    p.communicate()

def get_junctions(sams,args):
  prog = re.compile('^[MDNX=]$')
  outsams = {}
  z = 0
  outs = []
  for sam in sams:
    z+=1
    outsams[z] = sam
    v = [x for x in sam.value('cigar_array') if prog.match(x['op'])]
    juncs = [i for i in range(0,len(v)) if v[i]['op'] =='N' and v[i]['val'] >= args.minimum_intron_size]
    for i in juncs:
      coord1 = sum([x['val'] for x in v[0:i]]) + sam.value('pos')
      coord2 = coord1 + v[i]['val'] 
      b1 = Bed(sam.value('rname'),coord1-2,coord1-1)
      b2 = Bed(sam.value('rname'),coord2-1,coord2)
      outs.append([Junction(b1,b2),z])
  return [outs,outsams]

if __name__=="__main__":
  main()
