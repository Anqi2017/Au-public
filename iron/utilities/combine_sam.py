#!/usr/bin/python
import argparse, os, sys, re
from subprocess import call, Popen, PIPE

def main():
  parser = argparse.ArgumentParser(description="Take multiple bam files and produce a single sorted bam output.")
  parser.add_argument('-o','--output',required=True,help="BAMFILE output name")
  parser.add_argument('--threads',type=int,default=1)
  parser.add_argument('--name',action='store_true',help="sort the BAM file by name")
  parser.add_argument('input',nargs='+',help="BAMFILE input file")
  args = parser.parse_args()
  for file in args.input:
    if not os.path.isfile(file):
      sys.stderr.write("ERROR: input file not found\n")
      sys.exit()
    m = re.search('(.*)\.[sb]am$',file)
    if not m:
      sys.stderr.write("ERROR: wrong input file format\n")
      sys.exit()
    sys.stderr.write("using: "+file+"\n")
  m2 = re.search('(.+)\.bam$',args.output)
  if not m2:
    sys.stderr.write("ERROR: output needs to be a .bam file\n")
    sys.exit()
  output_filebase = m2.group(1)
  #get the appropriate header first
  sq = {}
  for file in args.input:
    cmd = 'samtools view -H '+file
    if re.search('\.sam',file):
      cmd = 'samtools view -SH '+file
    with os.popen(cmd) as stream:
      for line in stream:
        m3 = re.match('@SQ\s+SN:(\S+)\s+LN:\d+',line)
        if m3:
          sq[m3.group(1)] = line
        elif re.match('@SQ',line):
          sys.stderr.write("Unsupported header SQ format: "+line+"\n")
          sys.exit()
  thread_option = ''
  if args.threads > 1: thread_option = ' -@ '+str(args.threads)+' '
  namestring = ''
  if args.name: namestring = '-n'
  p = Popen('samtools view -Sb - | samtools sort '+namestring+' '+thread_option+' - '+output_filebase,shell=True,stdin=PIPE)
  # for the first file use the header to make a new header
  cmd = 'samtools view -H '+args.input[0]
  if re.search('\.sam$',args.input[0]):
    cmd = 'samtools view -HS '+args.input[0]
  sq_cnt = 0
  with os.popen(cmd) as stream:
    for line in stream:
      if re.match('@SQ',line): 
        sq_cnt+=1
      else:
        p.stdin.write(line)
      if sq_cnt == 1: #only write the sequences once (the first time we encounter one)
        for name in sorted(sq.keys()):
          p.stdin.write(sq[name])
  for file in args.input:
    cmd = 'samtools view -S '+file
    if re.search('\.bam',file):
      cmd = 'samtools view '+file
    with os.popen(cmd) as stream:
      for line in stream:
        p.stdin.write(line)
  p.communicate()[0]
  #cmd = 'samtools view -Sb '+args.input+' | samtools sort - '+m.group(1)+'.sorted' 
  #sys.stderr.write(cmd+"\n")
  #call(cmd,shell=True)

if __name__=="__main__":
  main()
