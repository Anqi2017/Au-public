#!/usr/bin/python
import sys, argparse, re, gzip
from Bio.Format.GPD import GPDStream
from Bio.Range import BedStream, union_range_array
from Bio.Stream import MultiLocusStream

def main():
  parser = argparse.ArgumentParser(description="For every gpd entry (sorted) intersect it with bed depth (sorted)",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('gpd_input',help="GPD file")
  parser.add_argument('bed_depth_input',help="GPD file")
  parser.add_argument('-o','--output',help="output file")
  args = parser.parse_args()
  
  inf1 = None
  if re.search('\.gz$',args.gpd_input):
    inf1 = gzip.open(args.gpd_input)
  else:
    inf1 = open(args.gpd_input)
  inf2 = None
  if re.search('\.gz$',args.bed_depth_input):
    inf2 = gzip.open(args.bed_depth_input)
  else:
    inf2 = open(args.bed_depth_input)
  gs = GPDStream(inf1)
  bs = BedStream(inf2)
  of = sys.stdout
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  mls = MultiLocusStream([gs,bs])
  z = 0
  for ml in mls:
    z += 1
    #if z%1000 == 0:
    sys.stderr.write(ml.get_range_string()+"       \r")
    [gpds,beds] = ml.get_payload()
    if len(gpds) == 0: 
      continue
    if len(beds)==0:
      for gpd in gpds:
        of.write(gpd.get_gene_name()+"\t"+gpd.get_transcript_name()+"\t"+str(gpd.get_exon_count())+"\t"+str(gpd.get_length())+"\t0\t0\t0"+"\n")
      continue
    #break beds up by depth
    #depths = {}
    #for bed in beds:
    #  d = int(bed.get_payload())
    #  if d not in depths: depths[d] = []
    #  depths[d].append(bed)
    #for gpd in gpds:
    #  clen = 0
    #  tot = 0
    #  for d in depths:
    #    covs = []
    #    for ex in [x.get_range() for x in gpd.exons]:
    #      clen += sum([x.overlap_size(ex) for x in depths[d]])
    #      tot += clen*d
    #  of.write(gpd.get_gene_name()+"\t"+gpd.get_transcript_name()+"\t"+str(gpd.get_exon_count())+"\t"+str(gpd.get_length())+"\t"+str(clen)+"\t"+str(float(clen)/float(gpd.get_length()))+"\t"+str(float(tot)/float(gpd.get_length()))+"\n")
    for gpd in gpds:
      covs = []
      for ex in [x.get_range() for x in gpd.exons]:
        c = union_range_array(ex,beds,payload=2)
        covs += c
      clen = sum([x.length() for x in covs if int(x.get_payload())>0])
      tot =  sum([x.length()*int(x.get_payload()) for x in covs])
      of.write(gpd.get_gene_name()+"\t"+gpd.get_transcript_name()+"\t"+str(gpd.get_exon_count())+"\t"+str(gpd.get_length())+"\t"+str(clen)+"\t"+str(float(clen)/float(gpd.get_length()))+"\t"+str(float(tot)/float(gpd.get_length()))+"\n")
  sys.stderr.write("\n")
  of.close()
  inf1.close()
  inf2.close()

if __name__=="__main__":
  main()
