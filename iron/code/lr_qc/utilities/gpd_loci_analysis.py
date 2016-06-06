#!/usr/bin/python
import argparse, sys, os, re, gzip
from shutil import rmtree
from multiprocessing import cpu_count, Pool, Lock
from tempfile import mkdtemp, gettempdir
from Bio.Format.GPD import GPDStream
from Bio.Stream import LocusStream
from Bio.Structure import TranscriptLoci, TranscriptLociMergeRules

glock = Lock()
last_range = None

def main():
  #do our inputs
  args = do_inputs()

  # Setup inputs 
  inf = sys.stdin
  if args.input != '-':
    if m.search('\.gz$',args.input):
      inf = gzip.open(args.input)
    else:
      inf = open(args.input)
  of = sys.stdout
  # Setup outputs
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')


  mr = TranscriptLociMergeRules('is_any_overlap')
  mr.set_use_junctions(False)
  if args.threads > 1:
    p = Pool(processes=args.threads)
  results = []
  z = 0
  for locus in LocusStream(GPDStream(inf)):
    if args.threads <= 1:
      tls = Queue(do_locus(locus,mr,z,args,verbose=True))
      results.append(tls)
    else:
      tls = p.apply_async(do_locus,args=(locus,mr,z,args,False),callback=process_output)
      results.append(tls)
    z += len(locus.get_payload())
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\n")
  sys.stderr.write("Outputing results\n")
  if args.output_loci:
    if re.search('\.gz$',args.output_loci):
      ofl = gzip.open(args.output_loci,'w')
    else:
      ofl = open(args.output_loci,'w')
  lnum = 0
  for res in sorted([y for y in [r.get() for r in results] if y],key=lambda x: (x.chr,x.start,x.end)):
    rng =  res.get_range_string()
    rngout = res.copy()
    tls = res.get_payload()
    for tl in sorted(tls,key=lambda x: (x.get_range().chr,x.get_range().start,x.get_range().end)):
      lnum += 1
      txs = sorted(tl.get_transcripts(),key=lambda x: (x.get_range().chr,x.get_range().start,x.get_range().end))
      tlrng = [str(x) for x in tl.get_range().get_bed_array()]
      ofl.write("\t".join(tlrng)+"\t"+str(lnum)+"\t"+str(len(txs))+"\n")
      for tx in txs:
        cov = tx.get_payload()[1]
        of.write("\t".join(tlrng)+"\t"+str(lnum)+"\t"+str(len(txs))+"\t"+str(tx.get_payload()[0])+"\t"+str(z)+"\t"+tx.get_gene_name()+"\t"+str(cov['average_coverage'])+"\t"+str(cov['fraction_covered'])+"\t"+str(cov['mindepth'])+"\n")
  if args.output_loci:
    ofl.close()
  inf.close()
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  #if not args.specific_tempdir:
  #  rmtree(args.tempdir)

def process_output(curr_range):
  global last_range
  global glock
  if not curr_range: return None
  if len(curr_range.get_payload())==0: return None
  glock.acquire()
  if curr_range:
    if not last_range: 
      last_range = curr_range
      sys.stderr.write("Pos: "+curr_range.get_range_string()+"          \r")
    elif curr_range.cmp(last_range) == 1:
      last_range = curr_range
      sys.stderr.write("Pos: "+last_range.get_range_string()+"          \r")
  glock.release()
  return curr_range

def do_locus(locus,mr,curline,args,verbose=True):
    txl = TranscriptLoci()
    txl.set_merge_rules(mr)
    #if len(locus.get_payload()) == 0: return
    z = 0
    tot = len(locus.get_payload())
    for gpd in locus.get_payload():
      z += 1
      curline+=1
      gpd.set_payload([curline,None])
      if verbose: sys.stderr.write("adding transcript: "+str(z)+'/'+str(tot)+"           \r")
      if gpd.get_exon_count() < args.min_exon_count: continue
      txl.add_transcript(gpd)
    if len(txl.get_transcripts()) > 0:
      covs = txl.get_depth_per_transcript()
      remove_list = []
      for g in txl.get_transcripts():
        if g.get_id() not in covs: continue
        x = covs[g.get_id()]
        g.get_payload()[1] = x
        if x['average_coverage'] < args.min_depth:
          remove_list.append(g.get_id())
        elif args.min_depth > 1 and x['fraction_covered'] < args.min_coverage_at_depth:
          remove_list.append(g.get_id())
      for tx_id in remove_list:
        txl.remove_transcript(tx_id)
    if verbose: sys.stderr.write("\n")
    if verbose: 
      if txl.get_range():
        sys.stderr.write(txl.get_range().get_range_string()+"\n")
    #sys.stderr.write('partition locus'+"\n")
    tls = txl.partition_loci(verbose=verbose)
    curr_range = None
    for tl in tls:
      if not curr_range: curr_range = tl.get_range()
      curr_range = curr_range.merge(tl.get_range())
    if not curr_range: return None
    curr_range.set_payload(tls)
    return curr_range

# Let us handle an output with get command like apply_async 
class Queue:
  def __init__(self,val):
    self.val = val
  def get(self):
    return self.val

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Read SORTED genepred as input. Output stuff about loci.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--output_loci',help="Only describe the loci")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Run specific arguments
  parser.add_argument('--min_depth',type=float,default=1,help="Only consider reads with averge coverage this much or higher")
  parser.add_argument('--min_coverage_at_depth',type=float,default=0.8,help="Only consider reads covered at 'min_depth' at this fraction or greater.")
  parser.add_argument('--min_exon_count',type=float,default=1,help="Only construct loci from reads with this many or more exons")

  ## Temporary working directory step 1 of 3 - Definition
  #group = parser.add_mutually_exclusive_group()
  #group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  #group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  # Temporary working directory step 2 of 3 - Creation
  #setup_tempdir(args)
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
