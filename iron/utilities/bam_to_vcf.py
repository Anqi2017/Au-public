#!/usr/bin/python
import argparse, sys, multiprocessing, re, os
from subprocess import Popen, PIPE

gcounter = 0
glock = multiprocessing.Lock()

def main():
  parser = argparse.ArgumentParser(description="Takes a BAM file preferably one already filtered to be uniquely mapped reads.")
  parser.add_argument('input_fasta',help="FASTAFILE indexed")
  parser.add_argument('input_sorted_bam',help="BAMFILE sorted indexed")
  parser.add_argument('--threads',type=int,default=multiprocessing.cpu_count(),help="Number of threads defautl cpu_count")
  parser.add_argument('--include_multiply_mapped_reads',action='store_true',help="Include multiply mapped reads that are excluded by default.  Note that this feature is not complete as it is with the 256 sam filter.  it will only remove secondary alignments while still leaving the multiply mapped primary alignments.  To only use uniquely mapped reads you need to pre-filter on unique and start from that indexed bam.")
  parser.add_argument('--include_indels',action='store_true',help="By default only SNPs and only loci with multiple genotypes are output.  This will output indels.")
  parser.add_argument('--consensus',action='store_true',help="Use the original caller")
  args = parser.parse_args()
  #read the sam header  
  p = Popen(('samtools view -H '+args.input_sorted_bam).split(),stdout=PIPE)
  chromlens = {}
  for line in p.stdout:
    m = re.match('@SQ\s+SN:(\S+)\s+LN:(\d+)',line.rstrip())
    if not m: continue
    chromlens[m.group(1)] = int(m.group(2))
  #Lets break these up now
  z = 0
  itersize = 10000000
  for chrom in chromlens:
    for i in range(1,chromlens[chrom],itersize):
      z+=1
  global gtotal
  gtotal = z
  if args.threads > 1:
    p = multiprocessing.Pool(processes=args.threads)
  for chrom in chromlens:
    for i in range(1,chromlens[chrom],itersize):
      rstart = i
      rend = itersize+i-1
      if rend > chromlens[chrom]: rend = chromlens[chrom]
      if args.threads <= 1:
        v = get_region_vcf(args,chrom,rstart,rend)
        do_output(v)
      else:
        p.apply_async(get_region_vcf,args=(args,chrom,rstart,rend),callback=do_output)
  if args.threads > 1:
    p.close()
    p.join()
def do_output(lines):
  global glock
  global gcounter
  global gtotal
  glock.acquire()
  gcounter += 1
  perc = "{0:.3g}%".format(100*float(gcounter)/float(gtotal))
  sys.stderr.write(perc+"        \r")
  for line in lines:
    print line.rstrip()
  glock.release()

def get_region_vcf(args,chrom,rstart,rend):
  region_arg = '-r '+chrom+':'+str(rstart)+'-'+str(rend)
  nullout = open(os.devnull,'w')
  excl_arg = ' --excl-flags 256 '
  if args.include_multiply_mapped_reads:
    excl_arg = ''
  vcf_arg = ' --skip-variants indels '
  if args.include_indels:
    vcf_arg = ''
  caller_arg = ' -mv '
  if args.consensus:
    caller_arg = ' -cv '
  p = Popen('samtools mpileup '+excl_arg+' '+region_arg+' -uf '+args.input_fasta+' '+args.input_sorted_bam+' | bcftools call '+caller_arg+' -Ov '+vcf_arg,shell=True,stdout=PIPE,stderr=nullout)
  lines = []
  for line in p.stdout:
    if line[0]=='#': continue
    lines.append(line.rstrip())
  return lines

if __name__=="__main__":
  main()
