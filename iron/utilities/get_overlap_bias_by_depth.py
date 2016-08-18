#!/usr/bin/python
import argparse, sys, os, gzip
from shutil import rmtree
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir

from Bio.Range import Bed

of = None

def do_open(fname):
  if fname[-3:]=='.gz':
    return gzip.open(fname)
  return open(fname)

def main():
  #do our inputs
  args = do_inputs()
  global of
  of = sys.stdout
  if args.output: of = open(args.output,'w')

  smallA = 9**9
  bigA = 0
  inf = do_open(args.depth_A)
  for line in inf:
    num = int(line.rstrip().split("\t")[3])
    if num < smallA: smallA = num
    if num > bigA: bigA = num
  inf.close()
  smallB = 9**9
  bigB = 0
  inf = do_open(args.depth_B)
  for line in inf:
    num = int(line.rstrip().split("\t")[3])
    if num < smallB: smallB = num
    if num > bigB: bigB = num
  inf.close()

  if args.threads > 1:
    p = Pool(processes=args.threads)
  for i in reversed(range(smallA,bigA+1)):
    for j in range(smallB,i+1):
      if args.threads > 1:
        p.apply_async(get_overlap,args=(args.depth_A,args.depth_B,i,j),callback=do_output)
      else:
        ov = get_overlap(args.depth_A,args.depth_B,i,j)
        do_output(ov)
  if args.threads > 1:
    p.close()
    p.join()
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_output(ov):
  global of
  of.write("\t".join([str(x) for x in ov])+"\t"+str(float(1)-float(ov[4])/float(ov[2]))+"\n")


def get_overlap(fileA,fileB,min_A,min_B):
  infA = do_open(fileA)
  infB = do_open(fileB)
  bufA = read_next(infA,min_A)
  bufB = read_next(infB,min_B)
  tot = 0
  sizeA = 0
  sizeB = 0
  if bufA:
    sizeA = bufA.length()
  if bufB:
    sizeB = bufB.length()
  zA = 1
  zB = 1
  while True:
    #if (zA%10000 ==0 or zB%10000==0): sys.stderr.write(str(zA)+" "+str(zB)+"  \r")
    if not bufA or not bufB: break
    c = bufA.cmp(bufB)
    if c == 0:
      tot += bufA.overlap_size(bufB)
      saveA = bufA
      nA = bufA.subtract(bufB)
      if len(nA) > 0 and nA[-1].end == bufA.end:
        num = bufA.get_payload()
        bufA = Bed(nA[-1].chr,nA[-1].start-1,nA[-1].end)
        bufA.set_payload(num)
      else:
        bufA = read_next(infA,min_A)
        if bufA:
          sizeA += bufA.length()
        zA+=1

      nB = bufB.subtract(saveA)
      if len(nB) > 0 and nB[-1].end == bufB.end:
        num = bufB.get_payload()
        bufB = Bed(nB[-1].chr,nB[-1].start-1,nB[-1].end)
        bufB.set_payload(num)
      else:
        bufB = read_next(infB,min_B)
        if bufB:
          sizeB += bufB.length()
        zB+=1

    elif c == -1:
      bufA = read_next(infA,min_A)
      if bufA:
        sizeA += bufA.length()
      zA += 1
    else:
      bufB = read_next(infB,min_B)
      if bufB:
        sizeB += bufB.length()
      zB += 1
  #sys.stderr.write("\n")
  if bufA:
    while True:
      bufA = read_next(infA,min_A)
      if bufA: sizeA += bufA.length()
      else: break
  if bufB:
    while True:
      bufB = read_next(infB,min_B)
      if bufB: sizeB += bufB.length()
      else: break
  infA.close()
  infB.close()
  return [min_A,min_B,sizeA,sizeB,tot]
  
def read_next(inf,strata):
  while True:
    v = inf.readline()
    if not v: return False
    num = int(v.rstrip().split("\t")[3])
    if num < strata: continue
    arr = v.split("\t")
    res = Bed(arr[0],int(arr[1]),int(arr[2]))
    res.set_payload(int(num))
    return res

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('depth_A',help="INPUT stratified bed depth file A")
  parser.add_argument('depth_B',help="INPUT stratified bed depth file B")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

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
