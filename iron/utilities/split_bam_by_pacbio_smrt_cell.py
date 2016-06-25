#!/usr/bin/python
import sys, argparse, re, os, gzip
from Bio.Format.Sam import BAMFile, sort_header
from subprocess import Popen, PIPE
from multiprocessing import Lock


def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Input bam file")
  parser.add_argument('-o','--output',required=True,help="output directory")
  parser.add_argument('--suffix',help="string to add before .bam: smrtcellname.XXXXXX.sorted.bam")
  args = parser.parse_args()
  args.output = args.output.rstrip('/')
  nameprog = re.compile('^([^\/]+)\/\d+/')
  if not os.path.exists(args.output):  os.makedirs(args.output)
  bf = BAMFile(args.input)
  sorted_header_text = sort_header(bf.header_text)
  fhs = {}
  z = 0
  for e in bf:
    z += 1
    if z%1000==0: sys.stderr.write(str(z)+" reads  "+str(len(fhs.keys()))+" cells      \r")
    m = nameprog.match(e.value('qname'))
    mol = '_nonpacbio'
    if m:
      mol = m.group(1)
    ln = e.get_line()
    if mol not in fhs:
      fname = args.output+'/'+mol
      if args.suffix:
        fname += '.'+args.suffix
      fname += '.gz'
      of = gzip.open(fname,'w')
      fhs[mol] = [of,fname]
      fhs[mol][0].write(sorted_header_text)
    fhs[mol][0].write(ln+"\n")
  sys.stderr.write("\n")
  z = 0
  for mol in fhs:
    z+=1
    fhs[mol][0].close()
    ofname = fhs[mol][1][:-2]+'sorted'
    inf = gzip.open(fhs[mol][1])
    cmd1 = 'samtools view -Sb -'
    cmd2 = 'samtools sort - '+ofname
    p2 = Popen(cmd2.split(),stdin=PIPE)
    p1 = Popen(cmd1.split(),stdin=inf,stdout=p2.stdin)
    p2.communicate()
    p1.communicate()
    inf.close()
    of.close()
    os.remove(fhs[mol][1])
    sys.stderr.write(str(z)+'/'+str(len(fhs.keys()))+" finished        \r")
  sys.stderr.write("\n")
if __name__=="__main__":
  main()
