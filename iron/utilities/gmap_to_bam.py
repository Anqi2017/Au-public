#!/usr/bin/python
import argparse, sys, multiprocessing, random, os, re, subprocess, gzip
from shutil import rmtree, copy
from tempfile import mkdtemp, gettempdir

def main():
  parser = argparse.ArgumentParser(description="Launch GMAP and run it in a temp directory until output is finished.")
  parser.add_argument('input_fastq',help="FASTA/FASTQ file or - for STDIN")
  parser.add_argument('output_bam')
  parser.add_argument('--gmap_index',required=True,help="Path to gmap index (directory)")
  parser.add_argument('--threads',type=int,default=multiprocessing.cpu_count())
  parser.add_argument('--max_paths',type=int,help="Maximum number of paths to show.")
  parser.add_argument('--max_intron_length',type=int,help="Maximum length of intron.")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  setup_tempdir(args)

  args.gmap_index = args.gmap_index.rstrip('/')
  m = re.match('^(.*)$',args.gmap_index)
  if m:
    gmapindexpath = m.group(1)
  m = re.search('([^\/]+)$',args.gmap_index)
  if m:
    gmapindexname = m.group(1)
  else:
    sys.stderr.write("problem reading gmap index "+args.gmap_index+"\n")
    sys.exit()

  maxpathpart = ''
  if args.max_paths:
    maxpathpart = ' -n '+str(args.max_paths)+' '

  maxintronpart = ''
  if args.max_intron_length:
    maxintronpart = ' -K '+str(args.max_intron_length)+' '

  #If its streamed lets just cache our input
  if args.input_fastq == '-' or args.input_fastq[-3:]=='.gz': 
    inf = sys.stdin
    if args.input_fastq[-3:]=='.gz': inf = gzip.open(args.input_fastq)
    of = open(args.tempdir+'/input.fastq','w')
    for line in inf:
      of.write(line)
    of.close()
    inf.close()
    args.input_fastq = args.tempdir+'/input.fastq'
    sys.stderr.write("Inputs prepared in: "+args.input_fastq+"\n")
  outtype = 'samse'
  gmap_cmd = 'gmap --ordered -D '+gmapindexpath+' -f '+outtype+' -d '+gmapindexname+' -t '+str(args.threads)+' '+maxpathpart+maxintronpart+' '+args.input_fastq
  sys.stderr.write("executing:\n"+gmap_cmd+"\n")

  allbam = args.tempdir+'/all.bam'
  of = open(allbam,'w')
  ofe = open(args.tempdir+'/log','w')
  cmd2 = "samtools view -Sb -"
  p2 = subprocess.Popen(cmd2.split(),stdin=subprocess.PIPE,stdout=of)
  p1 = subprocess.Popen(gmap_cmd.split(),stdout=p2.stdin,stderr=ofe)
  p1.communicate()
  p2.communicate()
  ofe.close()
  of.close()
  with open(args.tempdir+'/log') as inf:
    cnt = 0
    for line in inf:
      if re.match('No paths found',line): cnt+=1
      else:
        sys.stderr.write(line)
    sys.stderr.write("No paths found for "+str(cnt)+" queries\n")
  copy(args.tempdir+'/all.bam',args.output_bam)
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 


if __name__=='__main__':
  main()
