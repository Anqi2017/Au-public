#!/usr/bin/python
import argparse, sys, os, re, gzip
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import PIPE, Popen
from uuid import uuid4

def main():
  #do our inputs
  args = do_inputs()

  id = uuid4().hex[:6]

  sys.stderr.write("Traversing bam to discover names\n")
  cmd = 'samtools view '+args.bam_input
  if args.aligned: cmd += ' -F 4'
  names = []
  p = Popen(cmd.split(),stdout=PIPE)
  prog = re.compile('^([^\t]+)')
  z = 0
  tof = open(args.tempdir+'/nlist.txt.gz','w')
  cmd1 = 'sort -S'+str(args.mem)+'G -R -k1,1 --parallel='+str(args.threads)+' -T '+args.tempdir
  cmd2 = 'gzip'
  sys.stderr.write(cmd1+"\n")
  po2 = Popen(cmd2.split(),stdout=tof,stdin=PIPE)
  po1 = Popen(cmd1.split(),stdout=po2.stdin,stdin=PIPE)
  for line in p.stdout:
    z += 1
    if z%10000==0: sys.stderr.write("traversed "+str(z)+"   \r")
    m = prog.match(line)
    po1.stdin.write(m.group(1)+"\t"+str(z)+"\n")
    #names.append([m.group(1),z])
  sys.stderr.write("\n")
  po1.communicate()
  po2.communicate()
  p.communicate()
  tof.close()
  ## Now we can get a set of names
  sys.stderr.write("Traversing numbers to select\n")
  prev = None
  inf = gzip.open(args.tempdir+'/nlist.txt.gz')
  buffer = ''
  cnt = 0
  tof = open(args.tempdir+'/nlist2.txt.gz','w')
  cmd2 = 'gzip'
  p2 = Popen(cmd2.split(),stdin=PIPE,stdout=tof)
  cmd1 = 'sort -S'+str(args.mem)+'G -n --parallel='+str(args.threads)+' -T '+args.tempdir
  sys.stderr.write(cmd1+"\n")
  p1 = Popen(cmd1.split(),stdin=PIPE,stdout=p2.stdin)
  for line in inf:
    f = line.rstrip().split("\t")
    if f[0] != buffer:
      buffer = f[0]
      cnt += 1
      if cnt > args.count: break
    p1.stdin.write(f[1]+"\n")
  inf.close()
  p1.communicate()
  p2.communicate()
  tof.close()
  of = sys.stdout
  if args.output:
    of = open(args.output+'.'+str(cnt-1)+'.'+id+'.bam','w')
    cmd = "samtools view -Sb -"
  else: cmd = 'cat'
  p0 = Popen(cmd.split(),stdout=of,stdin=PIPE,close_fds=True)
  cmd = "samtools view -H "+args.bam_input
  p1 = Popen(cmd.split(),stdout=p0.stdin)
  p1.communicate()
  cmd = 'samtools view '+args.bam_input
  if args.aligned: cmd += ' -F 4'
  p = Popen(cmd.split(),stdout=PIPE)
  z = 0
  inf = gzip.open(args.tempdir+'/nlist2.txt.gz')
  curr = inf.readline()
  if curr: curr = int(curr.rstrip())
  for line in p.stdout:
    if not curr: continue
    z += 1
    if z%10000==0: sys.stderr.write("Final traversal of "+str(z)+" alignments    \r")
    if curr == z:
      p0.stdin.write(line) # write to the p0 pipe
      curr = inf.readline()
      if curr: curr = int(curr.rstrip())
  inf.close()
  p.communicate()
  sys.stderr.write("\n")
  p0.communicate()  
  of.close()


  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('bam_input',help="BAM file")
  parser.add_argument('count',type=int,help="number of read names to select")
  parser.add_argument('--aligned',action='store_true',help="only consider aligned reads")
  parser.add_argument('-o','--output',help="Output file basename, program will add .<size>.<RNG>.bam will use STDOUT to SAM if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="processor count")
  parser.add_argument('--mem',type=int,default=3,help="number of gigs of memory for buffering")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

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
