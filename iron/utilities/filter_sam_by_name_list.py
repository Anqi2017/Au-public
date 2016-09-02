#!/usr/bin/python
import sys, argparse, re
from subprocess import PIPE, Popen

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="BAM file or - for STDIN")
  parser.add_argument('name_list',help="Names to filter on")
  parser.add_argument('-o','--output',help="Output to bam file, leave blank for stdout")
  parser.add_argument('--inv',action='store_true',help="Get all names NOT in this list")
  args = parser.parse_args()
  
  names = set()
  with open(args.name_list) as inf:
    for line in inf:
      name = line.rstrip()
      names.add(name)
  inf = sys.stdin
  if args.input != '-':
    cmd = 'samtools view -h '+args.input
    p = Popen(cmd.split(),stdout=PIPE)
    inf = p.stdout
  of = sys.stdout
  if args.output:
    cmd = 'samtools view -Sb - -o '+args.output
    po = Popen(cmd.split(),stdin=PIPE,close_fds=True)
    of  = po.stdin
  in_header = True
  prog = re.compile('^([^\t]+)')
  for line in inf:
    if in_header:
      f = line.rstrip().split("\t")
      if len(f) > 10: 
        in_header = False
      else:
        of.write(line)
        continue  # move on since its a header
    # not header if we are still in
    m = prog.match(line)
    if args.inv and not m.group(1) in names:
      of.write(line)
    elif m.group(1) in names:
      of.write(line)

  if args.input != '-':
    p.communicate()
  if args.output:
    po.communicate()
if __name__=="__main__":
  main()
