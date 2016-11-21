#!/usr/bin/python
import os, sys, gzip, re, argparse
from subprocess import PIPE, Popen
def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="BAM file (alignqc indexed)")
  parser.add_argument('-o','--output',help="bam output or - for sam stdout")
  args = parser.parse_args()

  if not os.path.exists(args.input+'.bgi'):
    sys.stderr.write("ERROR you need an alignqc index\n")
    sys.exit()
  
  fname = args.input+'.bgi'
  inf = gzip.open(fname)
  bases = {}
  z = 0
  for line in inf:
    z += 1
    if z%1000 == 0: sys.stderr.write(str(z)+" paths assessed\r")
    f = line.rstrip().split("\t")
    m = re.match('(.+)_[^_]+$',f[0])
    m1 = re.match('(.+)_2D$',f[0])
    m2 = re.match('(.+)_tem$',f[0])
    m3 = re.match('(.+)_com$',f[0])
    if not m: sys.exit()
    name = f[0]
    base = m.group(1)
    if base not in bases:
      bases[base] = {}
      bases[base]['2D'] = {}
      bases[base]['2D']['l'] = -1
      bases[base]['2D']['name'] = ''
      bases[base]['tem'] = {}
      bases[base]['tem']['l'] = -1
      bases[base]['tem']['name'] = ''
      bases[base]['com'] = {}
      bases[base]['com']['l'] = -1
      bases[base]['com']['name'] = ''
    l = int(f[4])
    if m1:
      if l > bases[base]['2D']['l']:
        bases[base]['2D']['l'] = l
        bases[base]['2D']['name'] = name
    if m2:
      if l > bases[base]['tem']['l']:
        bases[base]['tem']['l'] = l
        bases[base]['tem']['name'] = name
    if m3:
      if l > bases[base]['com']['l']:
        bases[base]['com']['l'] = l
        bases[base]['com']['name'] = name
  sys.stderr.write("\n")  
  # now we want to 
  best_names = set()
  for base in bases:
    if bases[base]['2D']['l'] > 0:
      best_names.add(bases[base]['2D']['name'])
    elif bases[base]['tem']['l'] > 0:
      best_names.add(bases[base]['tem']['name'])
    elif bases[base]['com']['l'] > 0:
      best_names.add(bases[base]['com']['name'])
    elif bases[base]['2D']['l'] == 0:
      best_names.add(bases[base]['2D']['name'])
    elif bases[base]['tem']['l'] == 0:
      best_names.add(bases[base]['tem']['name'])
    elif bases[base]['com']['l'] == 0:
      best_names.add(bases[base]['com']['name'])
    else: 
      sys.stderr.write("ERROR. what?\n")
      sys.exit()

  of = sys.stdout
  if args.output:
    cmd = 'samtools view -Sb - -o '+args.output
    po = Popen(cmd.split(),stdin=PIPE)
    of = po.stdin
  cmd = 'samtools view -H '+args.input
  p = Popen(cmd.split(),stdout=of)
  p.communicate()

  cmd = 'samtools view '+args.input
  p = Popen(cmd.split(),stdout=PIPE)
  for line in p.stdout:
    f = line.rstrip().split("\t")
    if f[0] in best_names:
      of.write(line)
  p.communicate()

  if args.output:
    po.communicate()
  else:
    of.close()

if __name__=="__main__":
  main()
