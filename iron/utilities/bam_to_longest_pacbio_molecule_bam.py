#!/usr/bin/python
import sys, argparse, os, re, gzip

from subprocess import PIPE, Popen


def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('-o','--output',help="bam (alignqc indexed) file, or stdout sam if not set")
  args = parser.parse_args()


  if not os.path.exists(args.input+'.bgi'):
    sys.stderr.write("ERROR. Please build an alignqc index for this file\n")
    sys.exit()
  inf = gzip.open(args.input+'.bgi')
  r1 = re.compile('([^\/]+\/\d+\/\S+)')
  r2 = re.compile('([^\/]+\/\d+)')
  bases = {}
  for line in inf:
    m1 = r1.match(line)
    if not m1:
      sys.stderr.write("ERROR not a pacbio read\n"+line+"\n")
      sys.exit()
    name = m1.group(1)
    m2 = r2.match(name)
    base = m2.group(1)
    if base not in bases:
      bases[base] = {}
      bases[base]['name'] = ''
      bases[base]['l'] = -1
    f = line.rstrip().split("\t")
    l = int(f[4])
    if l > bases[base]['l']:
      bases[base]['l'] = l
      bases[base]['name'] = name
  inf.close()
  best_names = set()
  for base in bases: best_names.add(bases[base]['name'])

  of = sys.stdout
  if args.output:
    cmd = 'samtools view -Sb - -o '+args.output
    po = Popen(cmd.split(),stdin=PIPE)
    of = po.stdin

  cmd = 'samtools view -H '+args.input
  p1 = Popen(cmd.split(),stdout=of)
  p1.communicate()
  cmd = 'samtools view '+args.input
  p2 = Popen(cmd.split(),stdout=PIPE)
  for line in p2.stdout:
    m = r1.match(line)
    if not m:
      sys.stderr.write("ERROR not a pacbio name\n"+line+"\n")
      sys.exit()
    if m.group(1) in best_names:
      of.write(line)
  p2.communicate()

  if args.output:
    po.communicate()
  else:
    of.close()

  


if __name__=="__main__":
  main()
