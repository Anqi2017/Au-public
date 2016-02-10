#!/usr/bin/python
import argparse, sys
from GenePredBasics import GenePredEntry
from RangeBasics import Bed, sort_ranges, merge_ranges, pad_ranges, subtract_ranges
from subprocess import Popen, PIPE
from PoissonBasics import probability_threshold

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('gpd_input')
  parser.add_argument('bam_input')
  parser.add_argument('--intergenic_buffer',default=10000,type=int)
  parser.add_argument('--window_size',default=10000,type=int)
  parser.add_argument('--bin_size',default=1000,type=int)
  parser.add_argument('--use_off_regions',action='store_true',help="Use a region even if there is no reads mapped to it.")
  parser.add_argument('--get_exons',action='store_true')
  args = parser.parse_args()
  chr_beds = {}
  gene_beds = []
  exon_beds = []
  sys.stderr.write("Reading genepred file\n")
  asum = 0
  atot = 0
  with open(args.gpd_input) as inf:
    for line in inf:
      g = GenePredEntry(line)
      asum += g.length()
      atot += 1
      grng = g.get_bed()
      grng.direction = None
      if grng.chr not in chr_beds:
        chr_beds[grng.chr] = grng.copy()
      chr_beds[grng.chr] = chr_beds[grng.chr].merge(grng)
      gene_beds.append(grng)
      for i in range(0,g.get_exon_count()):
        erng = Bed(g.value('chrom'),g.value('exonStarts')[i],g.value('exonEnds')[i])
        exon_beds.append(erng)
  avglen = float(asum)/float(atot)
  sys.stderr.write("Sorting gene bed\n")
  gene_beds = sort_ranges(gene_beds)
  gene_beds = merge_ranges(gene_beds,already_sorted=True)
  sys.stderr.write("Sorting chromosome beds\n")
  chr_beds = sort_ranges([chr_beds[x] for x in chr_beds.keys()])
  sys.stderr.write("Sorting exon beds\n")
  exon_beds = sort_ranges(exon_beds)
  sys.stderr.write("Get padded genes\n")
  padded_gene_beds = pad_ranges(gene_beds,args.intergenic_buffer,chr_beds)
  padded_gene_beds = merge_ranges(padded_gene_beds,already_sorted=True)
  sys.stderr.write("Get intergenic regions\n")
  intergenic_beds = subtract_ranges(chr_beds,padded_gene_beds,already_sorted=True)
  intergenic_beds = merge_ranges(intergenic_beds,already_sorted=True)
  intergenic_beds = window_break(intergenic_beds,args.window_size)
  #for i in intergenic_beds: print i.get_range_string()
  sys.stderr.write("Get merged exons\n")
  exon_beds = merge_ranges(exon_beds)
  sys.stderr.write("Get introns\n")
  intron_beds = subtract_ranges(gene_beds,exon_beds,already_sorted=True)  
  intron_beds = merge_ranges(intron_beds,already_sorted=True)
  intron_beds = window_break(intron_beds,args.window_size)
  sys.stderr.write("Going through short reads\n")
  cmd = "sam_to_bed_depth.py "+args.bam_input
  p = Popen(cmd.split(),stdout=PIPE)
  for x in intron_beds: x.set_payload([]) # payloads are read depths
  for x in intergenic_beds: x.set_payload([]) # payloads are read depths
  for x in exon_beds: x.set_payload([]) # payloads are read depths
  introndepth = []
  intergenicdepth = []
  exondepth = []
  pseudoreadcount = 0
  if not args.get_exons: exon_beds = []
  section_count = 0
  while True:
    section_count += 1
    line = p.stdout.readline()
    if not line: break
    f = line.split("\t")
    depth = int(f[3])
    curr = Bed(f[0],int(f[1]),int(f[2]))
    if section_count %100==0: sys.stderr.write(curr.get_range_string()+"          \r")
    pseudoreadcount += depth
    if len(exon_beds) > 0:
      while curr.cmp(exon_beds[0]) > 0 and len(exon_beds) > 0: # we've passed the region
        v = exon_beds.pop(0)
        if len(v.get_payload()) == 0 and not args.use_off_regions: continue
        av = average(v)
        exondepth.append(av)
        #print str(av)+" exonic "+v.get_range_string()
      c = curr.cmp(exon_beds[0])
      if c == 0: # overlaps with intron
        size = curr.overlap_size(exon_beds[0])
        for i in range(0,size): exon_beds[0].get_payload().append(depth)
    if len(intron_beds) > 0:
      while curr.cmp(intron_beds[0]) > 0 and len(intron_beds) > 0: # we've passed the region
        v = intron_beds.pop(0)
        if len(v.get_payload()) == 0 and not args.use_off_regions: continue
        av = average(v)
        introndepth.append(av)
        #print str(av)+" intronic "+v.get_range_string()
      c = curr.cmp(intron_beds[0])
      if c == 0: # overlaps with intron
        size = curr.overlap_size(intron_beds[0])
        for i in range(0,size): intron_beds[0].get_payload().append(depth)
    if len(intergenic_beds) > 0:
      while curr.cmp(intergenic_beds[0]) > 0 and len(intergenic_beds) > 0: # we've passed the region
        v = intergenic_beds.pop(0)
        if len(v.get_payload()) == 0 and not args.use_off_regions: continue
        av = average(v)
        intergenicdepth.append(av)
        display(curr,introndepth,intergenicdepth,pseudoreadcount,avglen)
        #print str(av)+" intergenic "+v.get_range_string()
      c = curr.cmp(intergenic_beds[0])
      if c == 0: # overlaps with intron
        size = curr.overlap_size(intergenic_beds[0])
        for i in range(0,size): intergenic_beds[0].get_payload().append(depth)
      #if c > 0: # we passed the intron
      #  v = intergenic_beds.pop(0)
      #  av = average(v)
      #  intergenicdepth.append(av)
      #  print str(av)+" intergenic "+v.get_range_string()
  if args.use_off_regions:
    for x in exon_beds: introndepth.append(average(x.get_payload()))
    for x in intron_beds: introndepth.append(average(x.get_payload()))
    for x in intergenic_beds: intergenicdepth.append(average(x.get_payload()))
  p.communicate()

def display(cbed,intron,intergenic,pseudoreadcount,avglen):
  #e = sorted(exon)
  io = sorted(intron)
  ig = sorted(intergenic)
  if len(io) == 0 or len(ig) == 0: return
  print cbed.get_range_string()
  print pseudoreadcount
  #print str(len(e))+"\t"+str(e[int(len(e)*0.05)]) + "\t" + str(e[int(len(e)*0.5)])+"\t"+ str(e[int(len(e)*0.95)])
  #print 'lambd95'+"\t"+str(probability_threshold(e[int(len(e)*0.05)],0.95))+"\t"+str(probability_threshold(e[int(len(e)*0.5)],0.95))+"\t"+str(probability_threshold(e[int(len(e)*0.95)],0.95))
  print str(len(io))+"\t"+str(io[int(len(io)*0.05)]) + "\t" + str(io[int(len(io)*0.5)])+"\t"+ str(io[int(len(io)*0.95)])
  #print 'lambd95'+"\t"+str(probability_threshold(io[int(len(io)*0.05)],0.95))+"\t"+str(probability_threshold(io[int(len(io)*0.5)],0.95))+"\t"+str(probability_threshold(io[int(len(io)*0.95)],0.95))
  print str(len(ig))+"\t"+str(ig[int(len(ig)*0.05)]) + "\t" + str(ig[int(len(ig)*0.5)])+"\t"+ str(ig[int(len(ig)*0.95)])
  #print 'lambd95'+"\t"+str(probability_threshold(ig[int(len(ig)*0.05)],0.95))+"\t"+str(probability_threshold(ig[int(len(ig)*0.5)],0.95))+"\t"+str(probability_threshold(ig[int(len(ig)*0.95)],0.95))
  print '--------'

def median(x):
  if len(x) == 0: return None
  return sorted(x)[len(x)/2]

def window_break(inranges,window_size):
  outputs = []
  if len(inranges) == 0: return outputs
  for inrange in inranges:
    start = inrange.start
    while start+window_size < inrange.end:
      outputs.append(Bed(inrange.chr,start,start+window_size-1))
      start += window_size
  return outputs
  
def average(rng):
  l = rng.length()
  tot = sum(rng.get_payload())
  avg = float(tot)/float(l)
  return avg

if __name__=="__main__":
  main()
