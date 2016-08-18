#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

from Bio.Format.GPD import GPD
from Bio.Range import Bed
from Bio.Statistics import average

def main():
  #do our inputs
  args = do_inputs()

  sys.stderr.write("Reading reference genepred\n")
  ref = {}
  tx_strand = {}
  z = 0
  with open(args.reference_genepred) as inf:
    for line in inf:
      gpd = GPD(line)
      gname = gpd.get_gene_name()
      tname = gpd.get_transcript_name()
      tx_strand[tname] = gpd.get_strand()
      if gname not in ref: ref[gname] = []
      ref[gname].append(gpd)
      z += 1
  sys.stderr.write("Read "+str(len(ref.keys()))+" genes and "+str(z)+" transcripts\n")

  if args.maximum_isoforms > 0:
    sys.stderr.write("Removing genes with more than "+str(args.maximum_isoforms)+" isoforms.\n")
    for gname in ref.keys():
      if len(ref[gname]) > args.maximum_isoforms: del ref[gname]
    sys.stderr.write("Now have "+str(len(ref.keys()))+" genes and "+str(sum([len(ref[x]) for x in ref.keys()]))+" transcripts\n")

  sys.stderr.write("Filtering by length "+str(args.minimum_length)+" bp\n")
  for gname in ref.keys():
    passing = []
    for gpd in ref[gname]:
      if gpd.get_length() < args.minimum_length: continue
      passing.append(gpd)
    if len(passing) == 0: del ref[gname]
    else: ref[gname] = passing
  sys.stderr.write("Now have "+str(len(ref.keys()))+" genes and "+str(sum([len(ref[x]) for x in ref.keys()]))+" transcripts\n")
  
  sys.stderr.write("Converting gpd into exon bed\n")
  beds = []
  for gname in ref.keys():
    for gpd in ref[gname]:
      tname = gpd.get_transcript_name()
      for i in range(0,len(gpd.exons)):
        ex = gpd.exons[i]
        beds.append(ex.get_range().get_bed_array()+[gname,tname,i])
  with open(args.tempdir+'/gpd.bed','w') as of:
    for bed in sorted(beds,key=lambda x: (x[0],x[1],x[2],x[3],x[4],x[5])):
      of.write("\t".join([str(x) for x in bed])+"\n")
  sys.stderr.write("intersecting with bed depth\n")
  of = open(args.tempdir+'/intersect.bed','w')
  cmd = 'bedtools intersect -wo -a - -b '+args.tempdir+'/gpd.bed'
  p = Popen(cmd.split(),stdin=args.bed_depth,stdout=of)
  p.communicate()
  coverage = {}
  sys.stderr.write("Reading the intersection\n")
  with open(args.tempdir+'/intersect.bed') as inf:
    for line in inf:
        f = line.rstrip().split("\t")
        gname = f[7]
        tname = f[8]
        depth = int(f[3])
        bed1 = Bed(f[0],int(f[1]),int(f[2]))
        bed2 = Bed(f[4],int(f[5]),int(f[6]))
        bed = bed1.union(bed2)
        bed.set_payload(depth)
        if gname not in coverage:
          coverage[gname] = {}
        if tname not in coverage[gname]:
          coverage[gname][tname] = []
        coverage[gname][tname].append(bed)
  transcript_depths = {}
  for gname in coverage:
    for tname in coverage[gname]:
      ref_gpd = [x for x in ref[gname] if x.get_transcript_name()==tname][0]
      rlen = ref_gpd.get_length()
      bases_covered = sum([x.length() for x in coverage[gname][tname]])
      bases_area = sum([x.length()*x.get_payload() for x in coverage[gname][tname]])
      avg_depth = float(bases_area)/float(rlen)
      if avg_depth < args.minimum_average_depth: continue
      if bases_covered < args.minimum_length: continue
      #print gname
      #print tname
      #print rlen
      #print bases_covered
      #print bases_area
      total_positions = {}
      for ex in ref_gpd.exons:
        b = ex.get_range().get_bed_array()
        for i in range(b[1],b[2]):
          total_positions[i] = 0 # zero indexed
      for b in coverage[gname][tname]:
        depth = b.get_payload()
        barr = b.get_bed_array()
        for i in range(barr[1],barr[2]):
          total_positions[i] = depth
      transcript_depths[tname] = total_positions
  sys.stderr.write("have information needed to plot from "+str(len(transcript_depths.keys()))+" transcripts\n")
  outputs = []
  for tname in transcript_depths:
    depths = transcript_depths[tname]
    positions = sorted(depths.keys())
    tx_len = len(positions)
    bins = {}
    for i in range(0,tx_len):
      bin = int(100*float(i)/float(tx_len))
      if bin not in bins: bins[bin] = []
      bins[bin].append(depths[positions[i]])
    for bin in bins:
      bins[bin] = average(bins[bin])
    biggest = float(max(bins.values()))
    tx_array = [float(bins[x])/biggest for x in sorted(bins.keys())]
    if tx_strand[tname] == '-':
      tx_array.reverse()
    #outputs.append(tx_array)
    args.output.write(tname+"\t"+"\t".join([str(x) for x in tx_array])+"\n")
  #for i in range(0,100):
  #  args.output.write("\t".join([str(x[i]) for x in outputs])+"\n")
  
  args.output.close()

  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('bed_depth',help="INPUT BED DEPTH FILE or '-' for STDIN")
  parser.add_argument('reference_genepred',help="INPUT reference genepred file")
  parser.add_argument('--transcript_list',help="Limit outputs to this list of transcripts")
  parser.add_argument('--maximum_isoforms',type=int,default=1,help="Don't process reference transcripts with genes with more than this number of reference transcripts, set to 0 or less for including all isoforms")
  parser.add_argument('--minimum_length',type=int,default=1500,help="Don't process genepred entries shorter than this length.\n")
  parser.add_argument('--minimum_average_depth',type=float,default=4,help="require at least this average depth")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  # Setup inputs 
  if args.bed_depth == '-':
    args.bed_depth = sys.stdin
  else:
    args.bed_depth = open(args.bed_depth)
  # Setup outputs
  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout
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
