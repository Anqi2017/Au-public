#!/usr/bin/python
import argparse, re, os, subprocess, sys
import SamBasics

def main():
  parser = argparse.ArgumentParser(description="Get read counts from sam or bam.")
  parser.add_argument('input',help="FILENAME sam or bam")
  parser.add_argument('--add_report',action='store_true',help="make a new file where we replace sam or bam with a .mapped_count")
  args = parser.parse_args()
  if args.add_report:
    m = re.match('(.+)\.[bs]am',args.input)
    if not m:
      sys.stderr.write("bad inputfile type should be .bam or .sam\n")
      sys.exit()
    baseinput = m.group(1)
  samtag = ''
  if re.search('\.sam$',args.input): samtag = '-S'
  z = 0
  #se = open('/dev/stderr','w')
  p = subprocess.Popen('sort | uniq | wc -l',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
  with os.popen('samtools view '+samtag+' '+args.input) as inf:
    for line in inf:
      z += 1
      if z%100000 ==0: 
        sys.stderr.write(str(z)+" alignments processed\r")
      line = line.rstrip()
      d = SamBasics.sam_line_to_dictionary(line)
      if not SamBasics.check_flag(d['flag'],4):
        if SamBasics.check_flag(d['flag'],64):
          p.stdin.write(d['qname']+'.1'+"\n")
        elif SamBasics.check_flag(d['flag'],128):
          p.stdin.write(d['qname']+'.2'+"\n")
        else:
          sys.stderr.write("Unrecognized\n")
          sys.exit()
  sys.stderr.write("\n")
  aligned_reads = int(p.communicate()[0].rstrip())
  if args.add_report:
    of = open(baseinput+'.mapped_reads','w')
    of.write(str(aligned_reads)+"\n")
    return
  print aligned_reads
  #if args.aligned:

if __name__=="__main__":
  main()
