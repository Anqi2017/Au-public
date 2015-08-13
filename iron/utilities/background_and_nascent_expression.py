#!/usr/bin/python
import argparse, os, sys, re, subprocess
from random import randint
from GenePredBasics import line_to_entry as genepred_line_to_entry
from SequenceBasics import read_fasta_into_hash
from shutil import rmtree
import BedToolsBasics

def main():
  parser = argparse.ArgumentParser(description="Analyze nascent RNA from transcriptomes.")
  parser.add_argument('-i','--input',required=True,help="FILENAME of alignment, - for STDIN")
  parser.add_argument('--input_type',default='sam',choices=['sam','bam','psl','bed','gpd'])
  parser.add_argument('-t','--transcriptome',help="GENEPRED FILE reference transcriptome.")
  parser.add_argument('-g','--genome',help="Genome file.")
  parser.add_argument('--exon_padding',type=int,default=100)
  parser.add_argument('--locus_padding',type=int,default=10000)
  parser.add_argument('--intergenic_bin_size',type=int,default=10000)
  parser.add_argument('--intronic_bin_size',type=int,default=1000)
  parser.add_argument('--top_expressing_bin_cutoff',type=float,default=0.1,help="Remove results in the top fraction of intergenic bins. Rather consider these mislabeled geneic regions.")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--specific_tempdir',help="DIRECTORY the exact directory to make (if necessary) and use")
  group.add_argument('--tempdir',default='/tmp',help="DIRECTORY that a temporary directory can be made in.")
  args = parser.parse_args()
  tdir = setup_tempdir(args)
  sys.stderr.write("working in "+tdir+"\n")
  # Get all exons from the transcriptome
  bounds = transcriptome_to_exons(args.transcriptome,tdir)
  if not args.genome:
    #If we didn't specify a genome lets just use some quasibounds
    of1 = open(tdir+'/genome_bounds.bed','w')
    of2 = open(tdir+'/genome_bounds.lengths','w')
    for chr in bounds:
      of1.write(chr+"\t"+str(bounds[chr][0])+"\t"+str(bounds[chr][1])+"\n")
      of2.write(chr+"\t"+str(bounds[chr][1])+"\n")
    of1.close()
    of2.close()
  #Make fatter exons to distance the introns from starts sites
  cmd = "bedtools slop -b "+str(args.exon_padding)+" -i "+tdir+'/merged_exons.bed -g '+tdir+'/genome_bounds.lengths > '+tdir+'/merged_padded_exons.bed'
  subprocess.call(cmd,shell=True)
  #Make fatter loci to distance the loci from intergenic
  cmd = "bedtools slop -b "+str(args.locus_padding)+" -i "+tdir+'/merged_loci.bed -g '+tdir+'/genome_bounds.lengths > '+tdir+'/merged_padded_loci.bed'
  subprocess.call(cmd,shell=True)
  #Get introns only
  cmd = "bedtools subtract -a "+tdir+'/merged_loci.bed -b '+tdir+'/merged_padded_exons.bed > '+tdir+'/introns.bed'
  subprocess.call(cmd,shell=True)
  #Get intergenic only
  cmd = "bedtools subtract -a "+tdir+'/genome_bounds.bed -b '+tdir+'/merged_padded_loci.bed > '+tdir+'/intergenic.bed'
  subprocess.call(cmd,shell=True)
  break_into_bins(tdir+'/intergenic.bed',tdir+'/intergenic_bins.bed',args.intergenic_bin_size)

  #Overlap bam file with the intergenic bins
  cmd = 'bedtools intersect -abam '+args.input+' -b '+tdir+'/intergenic_bins.bed -wo -bed -split > '+tdir+'/reads_intergenic_bin_intersect.bed'
  subprocess.call(cmd,shell=True)
  #Get nonzero contents of bins
  bins = process_bins(tdir+'/reads_intergenic_bin_intersect.bed')
  lambda_intergenic = calculate_lambda(bins,args,args.intergenic_bin_size)
  # get the number of reads in the experiment
  cmd = 'cut -f 4 '+tdir+'/reads_intergenic_bin_intersect.bed | sort | uniq | wc -l > '+tdir+'/intergenic_bins_read_count.txt'
  subprocess.call(cmd,shell=True)
  readcount = 0
  with open(tdir+'/intergenic_bins_read_count.txt') as inf:
    readcount = int(inf.readline().rstrip())
  intergenic_rpk_distro = get_rpk_distribution(bins)
  intergenic_rpkm_distro = get_rpkm_distribution(bins,readcount)
  print "Intergenic results:"
  print str(readcount) + "\tintergenic reads"
  print str(lambda_intergenic)+"\tintergenic lambda cutting top fraction of "+str(args.top_expressing_bin_cutoff)
  # Now lets process intronic bins
  break_into_bins(tdir+'/introns.bed',tdir+'/intronic_bins.bed',args.intronic_bin_size)
  cmd = 'bedtools intersect -abam '+args.input+' -b '+tdir+'/intronic_bins.bed -wo -bed -split > '+tdir+'/reads_intronic_bin_intersect.bed'
  subprocess.call(cmd,shell=True)
  intronic_bins = process_bins(tdir+'/reads_intronic_bin_intersect.bed')
  # get the number of reads in the experiment
  cmd = 'cut -f 4 '+tdir+'/reads_intronic_bin_intersect.bed | sort | uniq | wc -l > '+tdir+'/intronic_bins_read_count.txt'
  subprocess.call(cmd,shell=True)
  intronic_readcount = 0
  with open(tdir+'/intronic_bins_read_count.txt') as inf:
    intronic_readcount = int(inf.readline().rstrip())
  print str(intronic_readcount) + "\tintronic reads"
  intronic_rpk_distro = get_rpk_distribution(intronic_bins)
  intronic_rpkm_distro = get_rpkm_distribution(intronic_bins,intronic_readcount)
  #print intronic_rpk_distro
  #print intronic_rpkm_distro
  print "percentile\tintergenic_rpk\tintergenic_rpkm\tintronic_rpkm"
  for i in range(0,100):
    print str(i)+"\t"+\
          str(intergenic_rpk_distro[i][0])+"\t"+\
          str(intergenic_rpkm_distro[i][0])+"\t"+\
          str(intronic_rpk_distro[i][0])+"\t"+\
          str(intronic_rpkm_distro[i][0])
  if not args.specific_tempdir:
    rmtree(tdir)

def get_rpk_distribution(bins):
  sizes = []
  for bin in bins:
    lval = 0
    for fnum in bins[bin]:
      lval += fnum
    sizes.append(lval)
  sizes.sort()
  return [[sizes[int(float(x)*0.01*float(len(sizes)))]*1000,x] for x in range(0,100)] 

def get_rpkm_distribution(bins,total_reads):
  sizes = []
  for bin in bins:
    lval = 0
    for fnum in bins[bin]:
      lval += fnum
    sizes.append(lval)
  sizes.sort()
  return [[get_rpkm(sizes[int(float(x)*0.01*float(len(sizes)))],total_reads),x] for x in range(0,100)] 

def get_rpkm(reads_in_gene,total_reads):
  return 1000000000*float(reads_in_gene)/(float(total_reads))

def calculate_lambda(bins,args,windows_size):
  sizes = []
  for bin in bins:
    lval = 0
    for fnum in bins[bin]:
      lval += fnum
    sizes.append(lval)
  sizes.sort()
  valid_sizes = sizes[:-1*int(len(sizes)*args.top_expressing_bin_cutoff)]
  lamb = 0
  total = 0
  for num in valid_sizes:
    total += 1
    lamb += num
  return windows_size*lamb/total

def calculate_direct_threshold(bins,args,thresh):
  sizes = []
  for bin in bins:
    lval = 0
    for fnum in bins[bin]:
      lval += fnum
    sizes.append(lval)
  sizes.sort()
  valid_sizes = sizes[:-1*int(len(sizes)*args.top_expressing_bin_cutoff)]
  ind = int(thresh*len(valid_sizes))
  if ind == len(valid_sizes): ind -= 1
  return valid_sizes[ind]

def process_bins(infile):
  bins = {}
  with open(infile) as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      locus = f[12] +"\t" + f[13] + "\t"+f[14]
      if locus not in bins:
        bins[locus] = []
      if float(f[15]) > 0:
        # store the fraction of the read that is overlapped divided by the length of the region
        bins[locus].append((float(f[15])/(float(f[2])-float(f[1])))/(float(f[14])-float(f[13])))
  return bins

def break_into_bins(infile,outfile,binsize):
  #if not os.path.exists(tdir+'/intergenic_bins'):
  #  os.makedirs(tdir+'/intergenic_bins')
  of = open(outfile,'w')
  with open(infile) as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      chr = f[0]
      start = int(f[1])
      finish = int(f[2])
      if finish-start < binsize: continue
      mystart = start
      while mystart+binsize < finish:
        of.write(chr+"\t"+str(mystart)+"\t"+str(mystart+binsize)+"\n")
        mystart += binsize
  of.close()

def transcriptome_to_exons(fname,tdir):
  of1 = open(tdir+'/all_exons.bed','w')
  of2 = open(tdir+'/all_loci.bed','w')
  bounds = {}
  with open(fname) as inf:
    for line in inf:
      if re.match('^#',line): continue
      e = genepred_line_to_entry(line)
      for i in range(0,len(e['exonStarts'])):
        if e['chrom'] not in bounds:
          bounds[e['chrom']] = [100000000000,0]
        if e['exonStarts'][i] < bounds[e['chrom']][0]:
          bounds[e['chrom']][0] = e['exonStarts'][i]
        if e['exonEnds'][i] > bounds[e['chrom']][1]:
          bounds[e['chrom']][1] = e['exonEnds'][i]
        of1.write(e['chrom']+"\t"+str(e['exonStarts'][i])+"\t"+str(e['exonEnds'][i])+"\n")
      of2.write(e['chrom']+"\t"+str(e['txStart'])+"\t"+str(e['txEnd'])+"\n")
  of1.close()
  of2.close()
  # Get the compressed exons
  cmd = "bedtools sort -i "+tdir+'/all_exons.bed > '+tdir+'/all_exons.sorted.bed'
  subprocess.call(cmd,shell=True)
  cmd = "bedtools merge -i "+tdir+'/all_exons.sorted.bed > '+tdir+'/merged_exons.bed'
  subprocess.call(cmd,shell=True)
  cmd = "bedtools sort -i "+tdir+'/all_loci.bed > '+tdir+'/all_loci.sorted.bed'
  subprocess.call(cmd,shell=True)
  cmd = "bedtools merge -i "+tdir+'/all_loci.sorted.bed > '+tdir+'/merged_loci.bed'
  subprocess.call(cmd,shell=True)
  return bounds
  
def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    return args.specific_tempdir.rstrip('/')
  dirname = args.tempdir.rstrip('/')+'/nas.'+str(randint(1,100000000))
  if not os.path.exists(dirname):
    os.makedirs(dirname)
  return dirname
  

if __name__=="__main__":
  main()
