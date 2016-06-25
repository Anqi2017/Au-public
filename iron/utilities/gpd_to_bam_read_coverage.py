#!/usr/bin/python
import sys, argparse
from Bio.Format.GPD import GPDStream
from Bio.Format.Sam import SamtoolsBAMStream, BAMFile
from Bio.Stream import MultiLocusStream
from Bio.Range import ranges_to_coverage, union_range_array

def main():
  parser = argparse.ArgumentParser(description="Intersect a bam with a gpd file to give bam coverage of each gpd entry",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('sorted_bam',help="sorted bam file")
  args = parser.parse_args()
  
  if args.input == '-':
    args.input = sys.stdin
  else: args.input = open(args.input)
  #bs = BAMFile(args.sorted_bam)
  bs = SamtoolsBAMStream(args.sorted_bam)
  gs = GPDStream(args.input)
  mls = MultiLocusStream([gs,bs])
  for ml in mls:
    [gpds,bams]=ml.get_payload()
    print ml
    print len(gpds)
    print len(bams)
    # easy job if there is no coverage
    #if len(gpds)==0: continue
    #if len(bams)== 0:
    #  for gpd in gpds:
    #    print gpd.get_gene_name()+"\t"+gpd.get_transcript_name()+"\t"+str(gpd.get_exon_count())+"\t"+str(gpd.get_length())+"\t0"
    #  continue
    #allrngs = []
    #for rng in [[y.get_range() for y in x.get_target_transcript(68).exons] for x in bams]:
    #  allrngs+=rng
    #cov = ranges_to_coverage(allrngs)
    #for gpd in gpds:
    #  tot = 0
    #  for exon in [x.get_range() for x in gpd.exons]:
    #    vals = union_range_array(exon,cov,payload=2)
    #    tot += sum([x.get_payload()*x.length() for x in vals])
    #  print gpd.get_gene_name()+"\t"+gpd.get_transcript_name()+"\t"+str(gpd.get_exon_count())+"\t"+str(gpd.get_length())+"\t"+str(float(tot)/float(gpd.get_length()))
if __name__=="__main__":
  main()
