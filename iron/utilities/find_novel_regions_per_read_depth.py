#!/usr/bin/python
import argparse, sys, subprocess, re, os
import GenePredBasics
from shutil import rmtree
from random import randint

# Pre: Take a sorted BAM file and find novel regions at various depths.
#      It also can take read annotations from annotate_psl_with_gpd.py
#      But if no read annotation is given, then it will find loci regardless
#      of annotation state
# Post: Report the reassembled locus in bed format for general location or 
#       genepred location for more details.
#       Report is in a bed fromat
#       <chromosome> <start base-0> <end base-1> <minimum depth> <exons in locus> <locus name>
# Modifies: Creates a temp directory and uses bedtools


def main():
  parser = argparse.ArgumentParser(description='report regions that lack annotations')
  parser.add_argument('--read_annotations',help="FILENAME either rawoutput or bestoutput from annotate_psl_with_gpd")
  parser.add_argument('--bam',help="FILENAME of sorted bam file",required=True)
  parser.add_argument('--tempdir',default="/tmp",help="DIRECTORY of where temporary files can be stored")
  parser.add_argument('--depth',type=int,help="INT Instead of checking many depths only check this depth")
  parser.add_argument('--minintron',default=68,type=int,help="INT minimum size of intron default 68")
  parser.add_argument('--maxintron',default=100000,type=int,help="INT maximum size of intron default 100000")
  parser.add_argument('--gpdoutput',help="FILENAME store the genepred file created")
  parser.add_argument('--output','-o',help="FILENAME bed format output")
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--full',action='store_true',help="Exclude reads with full matches, retaining only partial and novel matches.")
  group2.add_argument('--partial',action='store_true',help="Exclude reads with partial matches, retaining only novel reads DEFAULT.")
  args = parser.parse_args()

  depth = {}
  if not os.path.exists(args.tempdir):
    sys.stderr.write("could not find temporary directory path\n")
    return
  if not os.path.exists(args.tempdir.rstrip("/")+"/weirathe"):
    os.makedirs(args.tempdir.rstrip("/")+"/weirathe")
  tdir = args.tempdir.rstrip("/") + "/weirathe/weirathe.orphan"+str(randint(1,10000000))
  sys.stderr.write("Using temporary directory: "+tdir+"\n")
  if not os.path.exists(tdir):
    os.makedirs(tdir)

  # iterate though read annotations
  annotated_reads = set()
  if args.read_annotations:
    with open(args.read_annotations) as inf:
      for line in inf:
        line = line.rstrip()
        if re.match('^psl_entry_id\s',line): continue
        if re.match('^$',line): continue
        f = line.split("\t")

        if args.full: # we only want the full matches
          if f[9] != 'Full': continue
        annotated_reads.add(f[1])

  if args.bam:
    # Later we will want to have chromosome lengths
    cmd0 = "samtools view -H "+args.bam
    ps0 = subprocess.Popen(cmd0.split(),stdout=subprocess.PIPE)
    of0 = open(tdir+"/lengths.txt",'w')
    for line in ps0.stdout:
      line = line.rstrip()
      if re.match('^@SQ',line):
        m1 = re.search('\sSN:(\S+)',line)
        m2 = re.search('\sLN:(\S+)',line)
        if m1 and m2: 
          of0.write(m1.group(1)+"\t"+m2.group(1)+"\n")
    of0.close()
    ps0.communicate()
    # first filter our bam
    cmd1 = "samtools view -h "+args.bam
    ps1 = subprocess.Popen(cmd1.split(),stdout=subprocess.PIPE)
    cmd2 = "samtools view -Sb -o "+tdir+"/temp.bam"+" -"
    ps2 = subprocess.Popen(cmd2.split(),stdin=subprocess.PIPE)
    for line in ps1.stdout:
      f = line.rstrip().split("\t")
      if len(f) < 9:
        ps2.stdin.write(line)
      if f[0] not in annotated_reads:
        ps2.stdin.write(line)
    ps1.stdout.close()
    ps2.communicate()

    # Now sort the new bam file
    cmd3 = "samtools sort "+tdir+"/temp.bam"+" "+tdir+"/temp.sorted"
    subprocess.call(cmd3.split())
    # Now get the coverage information
    cmd4 = "bedtools genomecov -bg -split -ibam "+tdir+"/temp.sorted.bam"
    coverage_file = tdir+"/temp.bed"
    of4 = open(coverage_file,'w')
    subprocess.call(cmd4.split(),stdout=of4)
    of4.close()
    #find our maxdepth
    maxdepth = 0
    with open(coverage_file) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        cov = int(f[3])
        if cov > maxdepth: maxdepth = cov
    print maxdepth

    # for all our depths make a bed file to explore 
    fhs = {}
    depths = []
    d = 1 #starting depth
    while d < maxdepth:
      depths.append(d)
      d*=2
    depths.append(maxdepth)
    if args.depth: depths = [args.depth]
    sys.stderr.write(str(depths)+"\n")
    for i in depths:
      fhs[i] = open(tdir+"/depth."+str(i)+".bed",'w')
    with open(coverage_file) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        cov = int(f[3])
        for i in depths:
          if cov >= i:
            fhs[i].write(line)
          else:
            continue
    for i in fhs:
      fhs[i].close()

    #sort the bed files
    for i in depths:
      cmd5 = "bedtools sort -i "+tdir+"/depth."+str(i)+".bed"
      of5 = open(tdir+"/depth."+str(i)+".sorted.bed",'w')
      subprocess.call(cmd5.split(),stdout=of5)
      of5.close

    # for each of our depths get the merged bed
    z = 0
    if args.gpdoutput:
      ofgpd = open(args.gpdoutput,'w')
    ofout = sys.stdout
    if args.output:
      ofout = open(args.output,'w')
    for i in depths:
      #compress_depth(tdir,i,args.minintron)
      bfile = tdir + "/depth."+str(i)+".sorted.bed"
      gpd_entries = GenePredBasics.bed_to_genepred(args.minintron,args.maxintron,bfile)
      for e in gpd_entries:
        z+=1
        iter = e.entry['name']
        name = "depth-"+str(i)+"_"+str(iter)
        e.entry['gene_name'] = str(i)
        e.entry['name'] = name
        line = e.get_line()
        length = e.length()
        exons = e.get_exon_count()
        if args.gpdoutput:
          ofgpd.write(line+"\n")
        ofout.write(e.entry['chrom'] + "\t" + str(e.entry['txStart']) + "\t" + str(e.entry['txEnd']) + "\t" + str(i) + "\t" + str(exons) + "\t" + str(length) + "\t" + name +"\n")

  rmtree(tdir)        

if __name__=="__main__":
  main()
