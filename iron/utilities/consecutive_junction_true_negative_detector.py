#!/usr/bin/python
import argparse, sys, os, gzip
from shutil import rmtree
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir

from subprocess import Popen, PIPE

from Bio.Format.GPD import GPDStream, GPD
from Bio.Stream import MultiLocusStream

def main():
  #do our inputs
  args = do_inputs()

  inf = sys.stdin
  if args.input:
    if args.input[-3:]=='.gz': inf = gzip.open(args.input)
    else: inf = open(args.input)

  of = open(args.tempdir+'/input.gpd.gz','w')
  sys.stderr.write("sorting our input\n")
  input_cnt = sort_gpd(inf,of,args)
  of.close()
  inf.close()

  rinf = None
  if args.reference[-3:] == '.gz':
    rinf = gzip.open(args.reference)
  else:
    rinf = open(args.reference)
  sys.stderr.write("sorting our reference\n")
  rof = open(args.tempdir+'/ref.gpd.gz','w')
  sort_gpd(rinf,rof,args)
  rof.close()

  # Now we can traverse the ordered files by locus
  inf_input = gzip.open(args.tempdir+'/input.gpd.gz')
  inf_ref = gzip.open(args.tempdir+'/ref.gpd.gz')
  
  gsi = GPDStream(inf_input)
  gsr = GPDStream(inf_ref)
  mls = MultiLocusStream([gsi,gsr])
  z = 0
  y = 0
  output = []
  if args.threads > 1:
    p = Pool(args.threads)
  sys.stderr.write("processing "+str(input_cnt)+" inputs\n")
  for rng in mls:
    z+=1
    if z%10 == 0: 
      perc = int(100*float(y)/float(input_cnt+1))
      sys.stderr.write(rng.get_range_string()+" "+str(y)+" inputs "+str(perc)+"%                      \r")
    (input_entries,reference_entries) = rng.get_payload()
    if len(input_entries)==0: continue
    # Lets convert these back to lines to make the easier to pass through multiprocessing
    igpds = [x.get_gpd_line() for x in input_entries]
    rgpds = [x.get_gpd_line() for x in reference_entries]
    y += len(input_entries)
    if args.threads == 1:
      output.append(MiniQueue(process_locus(igpds,rgpds,args)))
    else:
      output.append(p.apply_async(process_locus,args=(igpds,rgpds,args,)))
  sys.stderr.write("\n")

  if args.threads > 1:
    p.close()
    p.join()

  of = sys.stdout
  if args.output: 
    if args.output[-3:] == '.gz':
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  tn_cnt = 0
  for out in output:
    outlines = out.get()
    for line in outlines:
      f = line.rstrip().split("\t")
      if int(f[3]) != 0: tn_cnt+=1
      of.write(line+"\n")
  of.close()
  perc = '?'
  if input_cnt > 0:
    perc = int(100*float(tn_cnt)/float(input_cnt))
  sys.stderr.write("Found "+str(tn_cnt)+" "+str(perc)+"% Unsupported Transcripts\n")
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

class MiniQueue:
  def __init__(self,val):
    self.val = [val]
  def get(self):
    return self.val.pop()

def process_locus(igpds,rgpds,args):
  input_entries = [GPD(x) for x in igpds]
  reference_entries = [GPD(x) for x in rgpds]
  outlines = []
  injun = get_consecutive_junctions(input_entries,args)
  refjun = get_consecutive_junctions(reference_entries,args)
  allrefjuncs = [] # consolidate reference junctions
  for refgpdset in refjun:
    (refgpd,refjuncs) = refgpdset
    for refjunc in refjuncs: allrefjuncs.append(refjunc) # append all reference junctions
  #sys.stderr.write("Now check the overlap\n")
  for ingpdset in injun:
    (ingpd,juncs) = ingpdset
    # one gpd at a time
    unsupported_pairs = junction_match(juncs,allrefjuncs,args)
    ostr = ''
    ostr += ingpd.get_gene_name()+"\t"
    ostr += ingpd.get_transcript_name()+"\t"
    ostr += str(len(juncs))+"\t"
    ostr += str(len(unsupported_pairs))+"\t"
    ostr += ";".join([x[0].get_string()+"~~"+x[1].get_string() for x in unsupported_pairs])
    outlines.append(ostr)
  return outlines

# is juncs1 in juncs2?
# return the unsupported pairs
def junction_match(juncs1,juncs2,args):
  unsupported = []
  for junc1 in juncs1:
    found_match = False
    for junc2 in juncs2:
      if junc1[0].overlaps(junc2[0],tolerance=args.junction_tolerance):
        if junc1[1].overlaps(junc2[1],tolerance=args.junction_tolerance):
          found_match = True
          break
    if not found_match:
      unsupported.append(junc1)
  return unsupported

def get_consecutive_junctions(entries,args):
  results = []
  for gpd in entries:
    jsets = []
    for i in range(0,len(gpd.junctions)-1):
      jsets.append([gpd.junctions[i],gpd.junctions[i+1]])
    results.append([gpd,jsets])
  return results

def sort_gpd(inf,of,args):
  sort_cmd = 'sort -k3,3 -k5,5n -k6,6n --parallel='+str(args.threads)
  gzip_cmd = 'gzip'
  p1 = Popen(gzip_cmd.split(),stdout=of,stdin=PIPE)
  p2 = Popen(sort_cmd.split(),stdout=p1.stdin,stdin=PIPE)
  input_cnt = 0
  for line in inf:
    f = line.rstrip().split("\t")
    if int(f[8])>2: 
      p2.stdin.write(line)
      input_cnt += 1
  p2.communicate()
  p1.communicate()
  return input_cnt

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Report transcripts not supported by junction pairs from a reference.\nOutput format is:\n<gene name> <transcript name> <total junction pair count> <unsupported junction pair count> <unsupported junction pairs>",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT GPD transcriptome to scan for negatives")
  parser.add_argument('-r','--reference',required=True,help="INPUT GPD consecutive junctions to scan")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  parser.add_argument('-j','--junction_tolerance',type=int,default=0,help="INT junction tolerance")

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
