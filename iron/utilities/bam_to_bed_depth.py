#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree
from multiprocessing import cpu_count, Lock, Pool
from tempfile import mkdtemp, gettempdir
from Bio.Format.Sam import BAMFile
from Bio.Stream import LocusStream
from Bio.Range import ranges_to_coverage, sort_genomic_ranges

current = 0
glock = Lock()
results = {}
of = None

def main():
  #do our inputs
  args = do_inputs()
  bf = BAMFile(args.input)
  if not args.all_alignments:
    if args.index_path:
      bf.read_index(args.index_path)
    else:
      bf.read_index(args.index_path)
  ls = LocusStream(bf)
  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout
  global of
  of = args.output
  z = 0
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for entries in ls:
    bedarray = []
    #print len(entries.get_payload())
    for e in entries.get_payload():
      if not args.all_alignments and not e.indexed_as_primary_alignment(): continue
      if not e.is_aligned(): continue
      tx = e.get_target_transcript(min_intron=args.minimum_intron_size)
      for exon in tx.exons:
        bedarray.append(exon.rng.copy())
        #print exon.rng.get_range_string()
    if len(bedarray) == 0: continue
    if args.threads > 1:
      p.apply_async(get_output,args=(bedarray,z,),callback=do_output)
    else:
      r = get_output(bedarray,z)
      do_output(r)
    z += 1
  if args.threads > 1:
    p.close()
    p.join()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)
  args.output.close()

def do_output(outputs):
  global glock
  global current
  global of
  global results
  glock.acquire()
  oline = outputs[0]
  z = outputs[1]
  results[z] = oline

  while current in results:
    prev = current
    of.write(results[current])
    del results[prev]
    current += 1
  glock.release()

def get_output(bedarray,z):
  sarray = sort_genomic_ranges(bedarray[:])
  covs = ranges_to_coverage(bedarray)
  olines = ''
  for c in covs:
    olines += c.chr+"\t"+str(c.start-1)+"\t"+str(c.end)+"\t"+str(c.get_payload())+"\n"
  return [olines,z]

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Convert a bam file into a bed file with depth.  By default only uses the primary alignment from a .bgi index file.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT BAM FILE")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--all_alignments',action='store_true',help="by default only do the primary alignments")
  parser.add_argument('--minimum_intron_size',default=68,type=int,help="any gaps smaller than this we close")
  parser.add_argument('--index_path',help="specity index location if its not in default place")
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
