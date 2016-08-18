#!/usr/bin/python
import argparse, sys, os, gzip
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE
def main(args):
  inf = sys.stdin
  if args.input:
    if args.input[-3:]=='.gz':
      inf = gzip.open(args.input)
    else: inf = open(args.input)

  extra = ''
  if args.extra:
    with open(args.extra) as inf:
      for line in inf:
        extra = line.rstrip()
        break
  
  #do our inputs
  procs = []
  for i in range(0,args.threads):
    cmd = 'sort -T '+args.tempdir+' -o '+args.tempdir+'/'+str(i)+'.frac'
    cmd += ' '+extra
    #sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdin=PIPE,close_fds=True)
    procs.append(p)
  z = 0
  for line in inf:
    z+=1
    rem = z%args.threads
    procs[rem].stdin.write(line)      
  inf.close()
  for p in procs:
    p.communicate()
  cmd = 'sort -m '+' '.join([args.tempdir+'/'+str(x)+'.frac' for x in range(0,args.threads)])
  cmd += ' '+extra
  #sys.stderr.write(cmd+"\n")
  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz':
      of = gzip.open(args.output,'w')
    else: of = open(args.output,'w')
  
  p = Popen(cmd.split(),stdout=of)
  p.communicate()
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  parser.add_argument('--extra',help="read extra args from file")
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


def external_cmd(cmd):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv


if __name__=="__main__":
  args = do_inputs()
  main(args)
