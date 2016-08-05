#!/usr/bin/python
import argparse, os, sys, re
from subprocess import call, Popen, PIPE
from Bio.Format.Sam import BAMFile

def main():
  parser = argparse.ArgumentParser(description="Take multiple bam files and produce a single sorted bam output.")
  parser.add_argument('-o','--output',required=True,help="BAMFILE output name")
  #parser.add_argument('--threads',type=int,default=1)
  #parser.add_argument('--sort',action='store_true',help="sort output")
  #parser.add_argument('--name',action='store_true',help="sort the BAM file by name")
  parser.add_argument('--numeric_names',action='store_true',help="order by integer in name")
  parser.add_argument('input',nargs='+',help="BAMFILE input file")
  args = parser.parse_args()
  names = args.input
  if args.numeric_names:
    names = sorted(args.input,key=lambda x: int(re.search('(\d+)[^\d]*$',x).group(1)))
  for file in names:
    if not os.path.isfile(file):
      sys.stderr.write("ERROR: input file not found\n")
      sys.exit()
    m = re.search('(.*)\.[sb]am$',file)
    if not m:
      sys.stderr.write("ERROR: wrong input file format\n"+file+"\n")
      sys.exit()
    sys.stderr.write("using: "+file+"\n")
  #m2 = re.search('(.+)\.bam$',args.output)
  #if not m2 and not args.output == '-':
  #  sys.stderr.write("ERROR: output needs to be a .bam file\n")
  #  sys.exit()
  #output_filebase = m2.group(1)
  #get the appropriate header first
  #thread_option = ''
  #if args.threads > 1: thread_option = ' --threads '+str(args.threads)+' '
  of = sys.stdout
  if args.output != '-':
    of = open(args.output,'w')
  cmd = 'samtools view -Sb - '
  p1 = Popen(cmd.split(),stdin=PIPE,stdout=of)
  cmd = 'samtools view -H '+args.input[0]
  p2 = Popen(cmd.split(),stdout=p1.stdin)
  p2.communicate()
  for file in names:
    cmd = 'samtools view '+file
    p3 = Popen(cmd.split(),stdout=p1.stdin)
    p3.communicate()
    sys.stderr.write('done '+file+"\n")
  p1.communicate()
  of.close()
  #cmd = 'samtools view -Sb '+args.input+' | samtools sort - '+m.group(1)+'.sorted' 
  #sys.stderr.write(cmd+"\n")
  #call(cmd,shell=True)

if __name__=="__main__":
  main()
