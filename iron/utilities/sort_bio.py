#!/usr/bin/python
import argparse, sys, os, re
from shutil import rmtree
#from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE
from Bio.Format.Sam import SamStream
from multiprocessing import cpu_count

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Sort these files by chromosome alphabetical, then start then end coordinate")
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--name',action='store_true',help="Sort by query name rather than location.  For GenePred this will default to gene name then the transcript name.")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--gpd',action='store_true')
  group2.add_argument('--bed',action='store_true')
  group2.add_argument('--psl',action='store_true')
  group2.add_argument('--bam',action='store_true',help="bam if file or sam if something else.")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()

  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args
# special case for the sam type
def do_sam(args):
  if args.input != '-':
    m = re.search('\.bam$',args.input)
    if not m:  
      sys.stderr.write("ERROR input expects bam unless piping to stdin.. then SAM with header\n")
      sys.exit()
  if not args.output:
    sys.stderr.write("ERROR sam sorts must output to a bam file\n")
    sys.exit()
  m = re.match('^(.+)\.bam$',args.output)
  if not m:
    sys.stderr.write("ERROR sam sorts must output to a bam file\n")
    sys.exit()
  cmdout = 'samtools sort - '+m.group(1)
  if args.threads:  cmdout += ' -@ '+str(args.threads)
  inf = None
  if args.input == '-':
    inf = sys.stdin
  else:
    cmd = 'samtools view -h '+args.input
    p = Popen(cmd.split(),stdout=PIPE,bufsize=1)
    inf = p.stdout
  s = SamStream(inf)
  header = s.header_text.rstrip().split("\n")
  split_stream = [header[i].split("\t") for i in range(0,len(header))] 
  sq_inds = [i for i in range(0,len(split_stream)) if split_stream[i][0]=='@SQ']
  nonsq_inds = [i for i in range(0,len(split_stream)) if split_stream[i][0]!='@SQ']
  top = [header[i] for i in nonsq_inds]
  chroms = sorted([split_stream[i] for i in sq_inds],key = lambda x: x[1][3:])
  cmd2 = 'samtools view -Sb -'
  pout = Popen(cmdout.split(),stdin=PIPE)
  p2 = Popen(cmd2.split(),stdin=PIPE,stdout=pout.stdin)
  for t in top:
    p2.stdin.write(t.rstrip()+"\n")
  for c in chroms:
    p2.stdin.write("\t".join(c).rstrip()+"\n")
  for sam in s:
    p2.stdin.write(sam.get_line().rstrip()+"\n")
  p2.communicate()
  pout.communicate()
  if args.input != '-':
    p.communicate()
  return

def main():
  #do our inputs
  args = do_inputs()
  # Temporary working directory step 3 of 3 - Cleanup
  #sys.stderr.write("working in: "+args.tempdir+"\n")
  if args.bam:
    do_sam(args)
    return
  cmd = "sort -T "+args.tempdir+'/'
  if args.psl:
    if args.name:
      cmd = "sort -k10,10 -T "+args.tempdir+'/'
    else:
      cmd = "sort -k14,14 -k15,15n -k16,16n -k9,9 -T "+args.tempdir+'/'
  if args.bed:
    cmd = "sort -k1,1 -k2,2n -T "+args.tempdir+'/'
  if args.gpd:
    if args.name:
      cmd = "sort -k1,1 -k2,2 -T "+args.tempdir+'/'
    else:
      cmd = "sort -k3,3 -k5,5n -k6,6n -k4,4 -T "+args.tempdir+'/'
  # Setup inputs 
  if args.input == '-':
    args.input = sys.stdin
  else:
    args.input = open(args.input)
  # Setup outputs
  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout
  p = Popen(cmd.split(),stdout=args.output,stdin=PIPE)
  for line in args.input:
    p.stdin.write(line)
  #p.stdin.close()
  #p.wait()
  p.communicate()
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

if __name__=="__main__":
  main()
