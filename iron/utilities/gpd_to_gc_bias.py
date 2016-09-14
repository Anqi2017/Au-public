#!/usr/bin/python
import argparse, sys, os, random, gzip
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir

from subprocess import PIPE, Popen

from Bio.Format.GPD import GPDStream
from Bio.Format.Fasta import FastaData
from Bio.Sequence import Seq

# The stratified best_X_covered option requires external calls and bedtools

def main():
  #do our inputs
  args = do_inputs()

  of = sys.stdout
  if args.output: of = open(args.output,'w')

  inf = sys.stdin
  if args.input != '-':
    if args.input[-3:] == '.gz':
      inf = gzip.open(args.input)
    else:
      inf = open(args.input)

  sys.stderr.write("reading in fasta\n")
  f = FastaData(open(args.reference).read())
  sh = GPDStream(inf)
  gc_bins = range(0,args.number_of_bins)
  bin_handles = []
  for i in range(0,args.number_of_bins):
    fname = args.tempdir+'/'+str(i)+'.bed.gz'
    cmd2 = 'bed_to_bed_depth.py - -o '+fname
    p2 = Popen(cmd2.split(),stdin=PIPE,close_fds=True)
    cmd1 = 'sort -k 1,1 -k2,2n -k3,3n -T '+args.tempdir
    p1 = Popen(cmd1.split(),stdin=PIPE,stdout=p2.stdin,close_fds=True)
    bin_handles.append([p1,p2,fname,i])

  if args.best_X_covered:
    sys.stderr.write("work out stratified data\n")
    cmd3 = 'bed_depth_to_stratified_coverage.py --minimum_coverage 10 --output_key '+args.tempdir+'/key'+' -r '+args.reference+' - -o '+args.tempdir+'/combo.bed.gz'
    pstrat3 = Popen(cmd3.split(),stdin=PIPE,close_fds=True)
    cmd2 = 'bed_to_bed_depth.py -'
    pstrat2 = Popen(cmd2.split(),stdin=PIPE,stdout=pstrat3.stdin,close_fds=True)
    cmd1 = 'sort -k 1,1 -k2,2n -k3,3n -T '+args.tempdir
    pstrat1 = Popen(cmd1.split(),stdin=PIPE,stdout=pstrat2.stdin,close_fds=True)

  num = 0
  for gpd in sh:
    num +=1
    if(num%1000==0): sys.stderr.write(str(num)+"     \r")
    results = []
    if args.minimum_sequence_length:
      if gpd.get_length() < args.minimum_sequence_length: continue
      seq = gpd.get_sequence(f).upper()
      seq_obj = Seq(seq)
      n_count = seq_obj.n_count()
      if len(seq) - n_count < args.min_non_N: continue
      gc = seq_obj.gc_content()
      gc_bin = int(args.number_of_bins*gc)
      if gc_bin == args.number_of_bins: gc_bin -=1;
      for exon in gpd.exons:
        bed_bin = ["\t".join([str(x) for x in  exon.rng.get_bed_array()]),gc_bin]
        results.append(bed_bin)
    elif args.fragment:
      seqlen = gpd.get_length()
      if seqlen < args.fragment: continue
      sfrags = int(float(seqlen)/float(args.fragment))
      sremain = seqlen % args.fragment
      offset = 0
      if random.random() < 0.5: offset = sremain
      #print '^^^'
      for i in range(0,sfrags):
        gsub = gpd.subset(i*args.fragment+offset,(i+1)*args.fragment+offset)
        seq = gsub.get_sequence(f).upper()
        seq_obj = Seq(seq)
        n_count = seq_obj.n_count()
        if len(seq)-n_count < args.min_non_N: continue
        gc = seq_obj.gc_content()
        gc_bin = int(args.number_of_bins*gc)
        if gc_bin == args.number_of_bins: gc_bin -=1;
        for exon in gsub.exons:
          bed_bin = ["\t".join([str(x) for x in  exon.rng.get_bed_array()]),gc_bin]
          results.append(bed_bin)
      
    for val in results:
      [bed,bin] = val
      bin_handles[bin][0].stdin.write(bed+"\n")
      if args.best_X_covered:
        pstrat1.stdin.write(bed+"\n")
        #if not gc: print len(gpd.get_sequence(f))
  
  sys.stderr.write("\n")
  for v in bin_handles:
    v[0].communicate()
    v[1].communicate()
  if args.best_X_covered:
    pstrat1.communicate()
    pstrat2.communicate()
    pstrat3.communicate()
    # If we want stratified data we should do it here
    sys.stderr.write("read the key\n")
    d = {}
    with open(args.tempdir+'/key') as inf:
      header = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        d[int(f[0])] = int(f[1])
    if args.best_X_covered not in d:
      sys.stderr.write("ERROR: the number of bases you specified is probably too big you didn't make the digit begin with 1 or 5 and restof the numbers be zero\n")
      sys.exit()
    num = d[args.best_X_covered]
    ninf = gzip.open(args.tempdir+'/combo.bed.gz')
    nof = gzip.open(args.tempdir+'/strat.bed.gz','w')
    for line in ninf:
      f = line.rstrip().split("\t")
      if int(f[3]) >= num:
        nof.write("\t".join(f[:-1])+"\n")
    nof.close()
    ninf.close()
    for i in range(0,len(bin_handles)):
      v = bin_handles[i]
      fname = v[2]
      fname2 = args.tempdir+'/'+str(v[3])+'.strata.bed.gz'
      gof = open(fname2,'w')
      cmd2 = 'gzip'
      p2 = Popen(cmd2.split(),stdout=gof,stdin=PIPE)
      cmd1 = 'bedtools intersect -a '+fname+' -b '+args.tempdir+'/strat.bed.gz'
      p1 = Popen(cmd1.split(),stdout=p2.stdin)
      p1.communicate()
      p2.communicate()
      gof.close()
      # lets just replace the name of the file that the final output will read from
      bin_handles[i][2]=fname2
  # Now we have bed depths for each bin
  for v in bin_handles:
    fname = v[2]
    #sys.stderr.write(fname+" ... prosessing\n")
    depths = {}
    bin = v[3]
    inf = gzip.open(fname)
    for line in inf:
      f = line.rstrip().split("\t")
      bases = int(f[2])-int(f[1])
      depth = int(f[3])
      if depth not in depths: depths[depth] = 0
      depths[depth] += bases
    inf.close()
    for depth in sorted(depths.keys()):
      of.write(str(bin)+"\t"+str(depth)+"\t"+str(depths[depth])+"\n")
  of.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Report the GC Bias of reads for overall data",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-r','--reference',required='True',help="Reference genome fasta")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  # Run specific options
  parser.add_argument('--min_non_N',type=int,default=50,help="at least this many non N bases in the sequence (or fragment if -f option)")
  parser.add_argument('-b','--number_of_bins',type=int,default=20,help="number of bins to separate gc content into")
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('-f','--fragment',type=int,help="Ignore shorter end fragments (randomly dropped)")
  group1.add_argument('--minimum_sequence_length',type=int,help="Ignore alignments shorter than this")

  parser.add_argument('--best_X_covered',type=int,help="only consider the best X covered bases.  leftmost digit must be a 1 or 5 and the rest zeros ie 1000 or 50000")

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
