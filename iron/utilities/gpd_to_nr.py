#!/usr/bin/python
import argparse, sys, os, gzip, random
from shutil import rmtree
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir

from Bio.Format.GPD import GPDStream, GPD
from Bio.Stream import LocusStream
from Bio.Structure import TranscriptGroup

def main(args):
  #do our inputs

  inf = sys.stdin
  if args.input != '-':
    if args.input[-3:] == '.gz':
      inf = gzip.open(args.input)
    else:
      inf = open(args.input)

  of = sys.stdout
  if args.output: 
    if args.output[-3:]=='.gz':
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  

  gs = GPDStream(inf)
  ls = LocusStream(gs)
  if args.threads > 1:
    p = Pool(processes=args.threads)
  results = []
  for locus_rng in ls:
    if args.threads == 1:
      sys.stderr.write(locus_rng.get_range_string()+"\n")
    else:
      sys.stderr.write(locus_rng.get_range_string()+"                 \r")
    gpds = locus_rng.get_payload()
    if args.threads > 1:
      new_gpds = p.apply_async(do_multi_round_locus,args=(gpds,args,))
      results.append(new_gpds)
    else:
      new_gpds = MiniQueue(do_multi_round_locus(gpds,args))
      results.append(new_gpds)
  if args.threads > 1:
    p.close()
    p.join()
  
  for result in results:
    new_gpds = result.get()
    for v in new_gpds:
      if not v['tx'].validate(): 
        sys.stderr.write("ERROR: invalid gpd entry\n")
        sys.stderr.write(v['tx'].get_fake_gpd_line()+"\n")
        sys.exit()
      fake_gpd = v['tx'].get_fake_gpd_line()
      #print v['tx'].get_gene_name()
      if args.gene_names: 
        f = fake_gpd.rstrip().split("\t")
        f[0] = v['tx'].get_gene_name()
        fake_gpd = "\t".join(f)
      of.write(fake_gpd+"\n")
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

class MiniQueue:
  def __init__(val):
    self.val = [val]
  def get():
    return self.val.pop(0)

def do_multi_round_locus(gpds,args):
    if args.threads == 1: sys.stderr.write("processing "+str(len(gpds))+" gpds\n")
    new_gpds = process_locus(gpds,args)
    if args.threads == 1: sys.stderr.write("merged to "+str(len(new_gpds))+" gpds\n")
    num_gpds = -1
    round = 1
    while num_gpds != len(new_gpds):
      round +=1
      num_gpds = len(new_gpds)
      buffer = []
      for v in new_gpds:
        #if v['evidence'] < args.minimum_support: continue
        for i in range(0,min(v['evidence'],max(args.minimum_support+1,args.minimum_junction_end_support+1))):
          nline = GPD(v['tx'].get_fake_gpd_line())
          # replace the gene name if we know it
          if not nline.validate():
            if args.threads == 1: sys.stderr.write("WARNING: 1. failed to make valid gpd. losing candidate\n")
            continue
          ngpd = GPD(nline.get_fake_gpd_line())
          if args.gene_names:
            ngpd.set_gene_name(v['tx'].get_gene_name())
          buffer.append(ngpd)
      gpds = buffer
      new_gpds = process_locus(gpds,args)
      if args.threads == 1: sys.stderr.write("round "+str(round)+" merged to "+str(len(new_gpds))+" gpds\n")
    return new_gpds

def process_locus(gpds,args):
  # get just gpds that are more than one exon
  gpds = [x for x in gpds if x.get_exon_count() > 1]
  if len(gpds) == 0: return []
  #sys.stderr.write('locus members: '+str(len(gpds))+"\n")
  # do longest members first as reference
  gpds = sorted(gpds, key=lambda x: x.get_exon_count(),reverse=True)

  if args.downsample and len(gpds) > args.downsample*2:
    #sys.stderr.write("warning: downsampling locus"+"\n")
    longest = gpds[0:args.downsample]
    remain = gpds[args.downsample:]
    random.shuffle(remain)
    gpds = longest+remain[0:args.downsample]
    #sys.stderr.write("Now working on: "+str(len(gpds))+"\n")
  # do a greedy combining of subsets
  gpd1 = gpds.pop(0)
  tx_group = TranscriptGroup()
  tx_group.add_transcript(gpd1,juntol=args.junction_tolerance)
  tx_groups = []
  tx_groups.append(tx_group)
  while len(gpds) > 0:
    if args.threads == 1: sys.stderr.write(str(len(gpds))+" remaining        \r")
    buffer = []
    for i in range(0,len(gpds)):
      added = tx_groups[-1].add_transcript(gpds[i],juntol=args.junction_tolerance,verbose=False)
      if added: continue
      buffer.append(gpds[i]) # save value not added
    for tx_group2 in tx_groups[0:-1]:  # make it non-greedy
      for tx2 in tx_group2.transcripts:
        added = tx_groups[-1].add_transcript(tx2,juntol=args.junction_tolerance,verbose=False)
    gpds = buffer
    if len(gpds) > 0:
      gpd1 = gpds.pop(0)
      tx_group = TranscriptGroup()
      tx_group.add_transcript(gpd1,juntol=args.junction_tolerance,verbose=False)
      tx_groups.append(tx_group)
  # now tx groups holds the merged down version
  if args.threads == 1: sys.stderr.write("\n")
  #sys.stderr.write('new membership: '+str(len(tx_groups))+"\n") 
  # Now we've finished this set
  z = 0
  buffer = []
  for tx_group in sorted(tx_groups,key=lambda x: len(x.transcripts),reverse=True):
    z += 1
    if len(tx_group.transcripts) < args.minimum_support: continue
    #sys.stderr.write("  "+str(z)+": "+str(len(tx_group.transcripts))+"\n")
    while True:
      # get rid of poorly supported ends according to lack of junctions
      if len(tx_group.junction_groups) <= 0: break
      if len(tx_group.junction_groups[0].evidence)<args.minimum_junction_end_support:
        tx_group.junction_groups.pop(0)
        #sys.stderr.write("prune left\n")
      elif len(tx_group.junction_groups[-1].evidence)<args.minimum_junction_end_support:
        tx_group.junction_groups.pop()
        #sys.stderr.write("prune right\n")
      else: break
    # See if we cut out everything
    if len(tx_group.junction_groups)==0: continue
    buffer.append(tx_group)
  tx_groups = buffer
  # now we have pruned groups
  #sys.stderr.write("after pruning: "+str(len(tx_groups))+"\n")
  results = []
  for tx_group in tx_groups:
    represent = tx_group.get_transcript()
    if args.gene_names: 
      represent.set_gene_name(tx_group.transcripts[0].get_gene_name())
    results.append({'tx':represent,'evidence':len(tx_group.transcripts)})
  return results

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Take a sorted gpd input and compress it down",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  # program commands
  parser.add_argument('--minimum_support',type=int,default=2,help="require at least this many reads supporting")
  parser.add_argument('--minimum_junction_end_support',type=int,default=2,help="require at least this many observations of an end junction")
  parser.add_argument('-j','--junction_tolerance',type=int,default=20,help="tolerance to consensus junction base to combine on")
  parser.add_argument('--downsample',type=int,default=250,help="keep this many of the longest, plus this many more at random")
  parser.add_argument('--gene_names',action='store_true',help="the genes are named appropriately already")

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
