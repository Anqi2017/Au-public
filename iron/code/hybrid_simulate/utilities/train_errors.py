#!/usr/bin/python
import argparse, sys, os, random, json, zlib, base64, gzip
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from Bio.Format.Sam import BAMFile
from Bio.Format.Fasta import FastaData
from collections import Counter
from Bio.Errors import ErrorProfileFactory

# Create an output 'error profile' object that contains
# Quality information
# Context error information
# General error information

def main(args):
  sys.stderr.write("Read reference fasta\n")
  fasta = FastaData(open(args.reference_fasta).read())
  sys.stderr.write("Read alignment file\n")
  bf = BAMFile(args.bam_input,reference=fasta)
  bf.read_index()
  total_qualities = []
  for j in range(0,100):
    total_qualities.append([])
  ef = ErrorProfileFactory()
  mincontext = 0
  alignments = 0
  for i in range(0,args.max_alignments):
    rname = random.choice(bf.index.get_names())
    coord = bf.index.get_longest_target_alignment_coords_by_name(rname)
    if not coord: continue
    bam = bf.fetch_by_coord(coord)
    qual = bam.value('qual')
    do_qualities(total_qualities,qual)
    if not bam.is_aligned(): continue
    alignments += 1
    ef.add_alignment(bam)
    if i%100 == 0:
      mincontext = ef.get_min_context_count('target')
      if mincontext:
        if mincontext >= args.min_context and alignments >= args.min_alignments: break
    sys.stderr.write(str(i+1)+" lines   "+str(alignments)+"/"+str(args.min_alignments)+" alignments   "+str(mincontext)+"/"+str(args.min_context)+" mincontext        \r")
  sys.stderr.write("\n")
  sys.stderr.write(str(mincontext)+" minimum contexts observed\n")
  target_context = ef.get_target_context_error_report()
  general_error_stats = ef.get_alignment_errors().get_stats()
  general_error_report = ef.get_alignment_errors().get_report()
  # convert report to table
  general_all = [x.split("\t") for x in general_error_report.rstrip().split("\n")]
  general_head = general_all[0]
  #print [y for y in general_all[1:]]
  general_data = [[y[0],y[1],int(y[2]),int(y[3])] for y in general_all[1:]]
  general_error_report = {'head':general_head,'data':general_data}
  quality_counts = []
  for vals in total_qualities:
    garr = []
    grp = {}
    for v in vals:
      if v[0] not in grp: grp[v[0]] = {}# check ordinal
      if v[1] not in grp[v[0]]: grp[v[0]][v[1]] = 0 # run length
      grp[v[0]][v[1]]+=1
    for ordval in sorted(grp.keys()):
      for runlen in sorted(grp[ordval].keys()):
        garr.append([ordval,runlen,grp[ordval][runlen]])
    quality_counts.append(garr)
  #Quailty counts now has 100 bins, each has an ordered array of
  # [ordinal_quality, run_length, observation_count]
  
  # Can prepare an output
  output = {}
  output['quality_counts'] = quality_counts
  output['context_error'] = target_context
  output['alignment_error'] = general_error_report
  output['error_stats'] = general_error_stats
  of = None
  if args.output[-3:]=='.gz':
    of = gzip.open(args.output,'w')
  else: of = open(args.output,'w')
  of.write(base64.b64encode(zlib.compress(json.dumps(output)))+"\n")
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_qualities(total_qualities,qual):
    qualities = []
    for j in range(0,100):
      qualities.append([])
    # break qual on homopolymers
    if not qual: return
    if len(qual) <= 1: return
    hp = [[qual[0]]]
    for i in range(1,len(qual)):
      if qual[i] ==hp[-1][0]: hp[-1]+=[qual[i]]
      else: hp += [[qual[i]]]
    ind = 0
    for vals in hp:
      frac = 100*float(ind)/float(len(qual))
      qualities[int(frac)].append([ord(vals[0]),len(vals)])
      ind += len(vals)
    prev = []
    for j in range(0,100):
      if len(qualities[j]) > 0: prev = qualities[j]
      else: qualities[j] = prev[:]
    for j in range(0,100):
      total_qualities[j]+=qualities[j]

def do_inputs():
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('bam_input',help="INPUT FILE")
  parser.add_argument('reference_fasta',help="Reference Fasta")
  parser.add_argument('-o','--output',required=True,help="OUTPUTFILE can be gzipped")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  parser.add_argument('--max_alignments',type=int,default=1000000,help="The absolute maximum number of alignments to try")
  parser.add_argument('--min_alignments',type=int,default=1000,help="Visit at least this many alignments")
  parser.add_argument('--min_context',type=int,default=10000,help="Stop after seeing this many of each context")
  
  
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

def external_cmd(cmd,version=None):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  #do our inputs
  args = do_inputs()
  main()
