#!/usr/bin/python
import sys, argparse, re
from GenePredBasics import GenePredEntry as GPD
from random import shuffle
from ArtificalReferenceSequenceBasics import ARS_conversion_string_factory as ACF, ARS
from RangeBasics import Bed
from SequenceBasics import read_fasta_into_hash

#pre: takes a genome fasta and a transcriptome
#     takes a number of fusions to output
#post: outputs all transcripts for which a simulated fusion lies in an intron


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('reference_genome')
  parser.add_argument('transcripts_genepred')
  parser.add_argument('--out_gpd',help="fusion genepred",required=True)
  parser.add_argument('--out_fasta',help="fusion fasta",required=True)
  parser.add_argument('--fusion_count',type=int,default=1000,help="Create this many fusions, max is number of genes/2.")
  args = parser.parse_args()
  ref = read_fasta_into_hash(args.reference_genome)
  of_gpd = open(args.out_gpd,'w')
  of_fasta = open(args.out_fasta,'w')
  genes = {}
  with open(args.transcripts_genepred) as inf:
    for line in inf:
      gpd = GPD(line.rstrip())
      if gpd.value('exonCount') <= 1: continue
      if gpd.value('gene_name') not in genes:
        genes[gpd.value('gene_name')] = []
      genes[gpd.value('gene_name')].append(gpd)
  gene_names = genes.keys()
  fusion_count = args.fusion_count
  shuffle(gene_names)
  pairs = []
  while True:
    if len(pairs) == fusion_count: break
    if len(gene_names) < 2: break
    pair = [gene_names[0],gene_names[1]]
    pairs.append(pair)
    gene_names.pop(0)
    gene_names.pop(0)
  for pair in pairs:
    [gpds,ars] = get_random_gpds_from_pair(pair,genes,ref)
    print ars.name
    of_fasta.write(ars.get_fasta())
    for gpd in gpds:
      of_gpd.write(gpd+"\n")
  of_gpd.close()
  of_fasta.close()

def get_random_gpds_from_pair(pair,genes,ref):
    #print 'gene 1 ('+pair[0]+'): '
    j1s = set()
    j1chrom = genes[pair[0]][0].value('chrom')
    j1starts = []
    j1ends = []
    j1strand = genes[pair[0]][0].value('strand')
    j2s = set()
    j2chrom = genes[pair[1]][0].value('chrom')
    j2starts = []
    j2ends = []
    j2strand = genes[pair[1]][0].value('strand')
    for gpd in genes[pair[0]]:
      if gpd.value('strand') != j1strand: continue
      if gpd.value('chrom') != j1chrom: continue
      j1starts.append(gpd.value('exonStarts')[0])
      j1ends.append(gpd.value('exonEnds')[-1])
      for j in gpd.calculate_junctions():
        j1s.add(j)
    #print 'gene 2 ('+pair[1]+'): '
    for gpd in genes[pair[1]]:
      if gpd.value('strand') != j2strand: continue
      if gpd.value('chrom') != j2chrom: continue
      j2starts.append(gpd.value('exonStarts')[0])
      j2ends.append(gpd.value('exonEnds')[-1])
      for j in gpd.calculate_junctions():
        j2s.add(j)
    j1shuf = list(j1s)
    shuffle(j1shuf)
    j2shuf = list(j2s)
    shuffle(j2shuf)
    #print j1shuf[0]
    #print j2shuf[0]
    if j1strand == '+':
      m = re.match('[^:]+:(\d+)',j1shuf[0])
      left = Bed(j1chrom,min(j1starts)-500,int(m.group(1))+500,j1strand)
      fsite1 = int(m.group(1))
    else:
      m = re.match('[^:]+:(\d+),[^:]+:(\d+)',j1shuf[0])
      left = Bed(j1chrom,int(m.group(2))-500,max(j1ends)+500,j1strand)
      fsite1 = int(m.group(2))
    if j2strand == '+':
      m = re.match('[^:]+:(\d+),[^:]+:(\d+)',j2shuf[0])
      right = Bed(j2chrom,int(m.group(2))-500,max(j2ends)+500,j2strand)
      fsite2 = int(m.group(2))
    else:
      m = re.match('[^:]+:(\d+),[^:]+:(\d+)',j2shuf[0])
      right = Bed(j2chrom,min(j2starts)-500,int(m.group(1))+500,j2strand)
      fsite2 = int(m.group(1))
    #print left.get_range_string()+' '+left.direction
    #print right.get_range_string()+' '+right.direction
    [leftcomp,rightcomp] = get_compatible_transcripts(genes[pair[0]],fsite1,genes[pair[1]],fsite2)
    #print fsite1
    #print fsite2
    acf = ACF()
    acf.add_bounds(left)
    acf.add_bounds(right)
    ln= leftcomp[0].value('gene_name')
    rn= rightcomp[0].value('gene_name')
    site_string = leftcomp[0].value('chrom')+":"+str(fsite1)+leftcomp[0].value('strand')+'/'+rightcomp[0].value('chrom')+":"+str(fsite2)+rightcomp[0].value('strand')
    ars = ARS(ref=ref,conversion_string=acf.get_conversion_string(),name=ln+","+rn+","+site_string)
    #print ars.conversion_string
    #print ars.name
    #print ars.get_ars_name()
    gpds = make_new_genepreds(leftcomp,fsite1,rightcomp,fsite2,ars)
    return [gpds,ars]

def make_new_genepreds(lefts,fsite1,rights,fsite2,ars):
  gpd_lines = []
  for l in lefts:
    for r in rights:
      #now lets traverse these putting them together
      coords_left = []
      real_left = []
      real_right = []
      detectable_left = None #the coordiante we can see
      detectable_right = None
      if l.value('strand')=='+':
        for i in range(0,l.get_exon_count()):
          if l.value('exonStarts')[i]+1 > fsite1: break #too far
          #print str(l.value('exonStarts')[i])+'-'+str(l.value('exonEnds')[i])
          start = ars.convert_genomic_to_ARS_coordinate(l.value('chrom'),l.value('exonStarts')[i]+1)
          end = ars.convert_genomic_to_ARS_coordinate(l.value('chrom'),l.value('exonEnds')[i])
          v= [start,end]
          coords_left.append(v)
          rv = [l.value('exonStarts')[i]+1,l.value('exonEnds')[i]]
          real_left.append(rv)
        detectable_left = real_left[-1][1]
      else:
        for i in range(0,l.get_exon_count()):
          if l.value('exonEnds')[i] < fsite1: continue #not far enough
          #print str(l.value('exonStarts')[i])+'-'+str(l.value('exonEnds')[i])
          start = ars.convert_genomic_to_ARS_coordinate(l.value('chrom'),l.value('exonStarts')[i]+1)
          end = ars.convert_genomic_to_ARS_coordinate(l.value('chrom'),l.value('exonEnds')[i])
          v = [start,end]
          coords_left.append(v)
          rv = [l.value('exonStarts')[i]+1,l.value('exonEnds')[i]]
          real_left.append(rv)
        detectable_left = real_left[0][0]
        coords_left.reverse()
        for x in coords_left: x.reverse()
      #print coords_left
      coords_right = []
      if r.value('strand')=='+':
        for i in range(0,r.get_exon_count()):
          if r.value('exonEnds')[i] < fsite2: continue #not far enough
          #print str(r.value('exonStarts')[i])+'-'+str(r.value('exonEnds')[i])+' '+r.value('strand')
          start = ars.convert_genomic_to_ARS_coordinate(r.value('chrom'),r.value('exonStarts')[i]+1)
          end = ars.convert_genomic_to_ARS_coordinate(r.value('chrom'),r.value('exonEnds')[i])
          v= [start,end]
          #print v
          coords_right.append(v)
          rv = [r.value('exonStarts')[i]+1,r.value('exonEnds')[i]]
          real_right.append(rv)
        detectable_right = real_right[0][0]
      else:
        for i in range(0,r.get_exon_count()):
          if r.value('exonStarts')[i]+1 > fsite2: break #too far
          #print str(r.value('exonStarts')[i])+'-'+str(r.value('exonEnds')[i])+' '+r.value('strand')
          start = ars.convert_genomic_to_ARS_coordinate(r.value('chrom'),r.value('exonStarts')[i]+1)
          end = ars.convert_genomic_to_ARS_coordinate(r.value('chrom'),r.value('exonEnds')[i])
          #print v
          v = [start,end]
          coords_right.append(v)
          rv = [r.value('exonStarts')[i]+1,r.value('exonEnds')[i]]
          real_right.append(rv)
        detectable_right = real_right[-1][1]
        coords_right.reverse()
        for x in coords_right: x.reverse()
      #print coords_right
      #print ars.get_conversion_string()
      exonstarts = []
      exonends = []
      for x in coords_left: exonstarts.append(x[0][0]-1)
      for x in coords_right: exonstarts.append(x[0][0]-1)
      for x in coords_left: exonends.append(x[1][0])
      for x in coords_right: exonends.append(x[1][0])
      site_string = l.value('chrom')+':'+str(fsite1)+l.value('strand')+'/'+r.value('chrom')+':'+str(fsite2)+r.value('strand')
      detectable_string = l.value('chrom')+':'+str(detectable_left)+l.value('strand')+'/'+r.value('chrom')+':'+str(detectable_right)+r.value('strand')
      ostr = ''
      ostr += ars.name+"\t"
      ostr += l.value('name')+','+r.value('name')+','+detectable_string+"\t"
      ostr += ars.get_ars_name()+"\t"
      ostr += '+'+"\t"
      ostr += str(exonstarts[0])+"\t"
      ostr += str(exonends[-1])+"\t"
      ostr += str(exonstarts[0])+"\t"
      ostr += str(exonends[-1])+"\t"
      ostr += str(len(exonstarts))+"\t"
      ostr += ','.join([str(x) for x in exonstarts])+','+"\t"
      ostr += ','.join([str(x) for x in exonends])+','
      gpd_lines.append(ostr)
  return gpd_lines
def get_compatible_transcripts(g1,f1,g2,f2):
  #do left
  site_string = g1[0].value('chrom')+':'+str(f1)+g1[0].value('strand')+'/'+g2[0].value('chrom')+':'+str(f2)+g2[0].value('strand')
  left_compatible = []
  right_compatible = []
  for gpd in g1:
    if gpd.value('strand') == '+':
      match = False
      for i in range(0,gpd.get_exon_count()-1):
        if gpd.value('exonEnds')[i] <= f1 and gpd.value('exonStarts')[i+1]+1 > f1:
          match = True
          break
      if match: left_compatible.append(gpd)
    if gpd.value('strand') == '-':
      match = False
      for i in range(0,gpd.get_exon_count()-1):
        if gpd.value('exonEnds')[i] < f1 and gpd.value('exonStarts')[i+1]+1 >= f1:
          match = True
          break
      if match: left_compatible.append(gpd)
  for gpd in g2:
    if gpd.value('strand') == '-':
      match = False
      for i in range(0,gpd.get_exon_count()-1):
        if gpd.value('exonEnds')[i] <= f2 and gpd.value('exonStarts')[i+1]+1 > f2:
          match = True
          break
      if match: right_compatible.append(gpd)
    if gpd.value('strand') == '+':
      match = False
      for i in range(0,gpd.get_exon_count()-1):
        if gpd.value('exonEnds')[i] < f2 and gpd.value('exonStarts')[i+1]+1 >= f2:
          match = True
          break
      if match: right_compatible.append(gpd)
  if len(left_compatible) == 0: 
    sys.stderr.write("ERROR: problem with left\n")
    sys.exit()
  if len(right_compatible) == 0: 
    sys.stderr.write("ERROR: problem with right\n")
    sys.stderr.write(site_string+"\n")
    for g in g2: sys.stderr.write(g.get_line()+"\n")
    sys.exit()
  return [left_compatible,right_compatible]          
if __name__=="__main__":
  main()
