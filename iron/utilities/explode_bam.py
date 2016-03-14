#!/usr/bin/python
import sys, argparse
from subprocess import Popen, PIPE
from SamBasics import SamStream
from multiprocessing import cpu_count, Pool
def main():
  parser = argparse.ArgumentParser(description="Break a bam into evenly sized chunks print the number of chunks",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN sam or directly name bamfile")
  parser.add_argument('output_base',help="output base name myout will go to myout.1.bam")
  parser.add_argument('-k',type=int,required=True,help="Number per chunk")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Number of threads")
  args = parser.parse_args()
  
  inf = None
  if args.input == '-':
    inf = sys.stdin
  else: 
    cmd = "samtools view -h "+args.input
    p = Popen(cmd.split(),stdout=PIPE)
    inf = p.stdout

  v = SamStream(inf)
  buffer = []
  i = 0
  if args.threads > 1:
    poo= Pool(processes=args.threads)
  while True:
    e = v.read_entry()
    if not e: break
    buffer.append(e)
    if len(buffer) >= args.k:
      i+=1
      if args.threads > 1:
        poo.apply_async(do_output,args=(buffer,v.header[:],i,args.output_base))
      else:
        do_output(buffer,v.header[:],i,args.output_base)
      buffer = []
  if len(buffer) > 0:
    i+=1
    if args.threads > 1:
      poo.apply_async(do_output,args=(buffer,v.header[:],i,args.output_base))
    else:
      do_output(buffer,v.header[:],i,args.output_base)
  if args.threads > 1:
    poo.close()
    poo.join()

  if args.input != '-':
    p.communicate()
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

if __name__=="__main__":
  main()
