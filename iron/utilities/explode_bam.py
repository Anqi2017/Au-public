#!/usr/bin/python
import sys, argparse, re, os
from subprocess import Popen, PIPE
from SamBasics import is_header
from multiprocessing import cpu_count, Pool

def main():
  parser = argparse.ArgumentParser(description="Break a bam into evenly sized chunks",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="name bam file")
  parser.add_argument('output_base',help="output base name myout will go to myout.1.bam")
  parser.add_argument('-k',type=int,required=True,help="Number per chunk")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Number of threads")
  parser.add_argument('--name',action='store_true',help="pre-sorted by query name keep queries together")
  parser.add_argument('-F',help="Add an input flag filter if you are reading from a bam file")
  args = parser.parse_args()
  
  # read header first
  header = []
  cmd = "samtools view -H "+args.input
  p = Popen(cmd.split(),stdout=PIPE,bufsize=1)
  inf = p.stdout
  for line in inf:
    header.append(line)
  p.communicate()

  inf = None
  cmd = "samtools view "+args.input
  if args.F:
    cmd += " -F "+args.F
  p = Popen(cmd.split(),stdout=PIPE,bufsize=1)
  rex = re.compile('^(\S+)')
  buffersize = args.k
  buffer = []
  prev_name = None
  i = 0
  poo = Pool(processes=max(1,args.threads-2))
  while True:
    line = p.stdout.readline()
    if not line: break
    if args.name:
      m = rex.match(line)
      if prev_name and m.group(1) != prev_name and len(buffer) >= buffersize:
        i+= 1 
        poo.apply_async(do_output,args=(buffer[:],header,i,args.output_base))
        buffer = []
      prev_name = m.group(1)
    else:
      if len(buffer) >= buffersize:
        i+=1
        poo.apply_async(do_output,args=(buffer[:],header,i,args.output_base))
        buffer = []
    buffer.append(line)
  # Deal with remainder
  if len(buffer) > 0: 
    i+=1
    poo.apply_async(do_output,args=(buffer[:],header,i,args.output_base))
    buffer = []
  poo.close()
  poo.join()
  print i

def do_output(buffer,header,i,output_base):
  of = open(output_base+'.'+str(i)+'.bam','w')
  cmd = 'samtools view - -Sb'
  p = Popen(cmd.split(),stdin=PIPE,stdout=of)
  for e in header:
    p.stdin.write(e)
  for e in buffer:
    p.stdin.write(e)
  p.communicate()
  of.close()
  return

if __name__=="__main__":
  main()
