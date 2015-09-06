#!/usr/bin/python
import sys

#######
#  Jason Weirather 20140317
#  This script will take genePred files exons and convert them to gtf files
#
#  Pre:  Input the filename of the genePred file
#  Post: Outputs a gtf format to standard output
#  Modifies:  Standard output
#######

if len(sys.argv) < 2:
  print 'gpd_to_gtf.py <Input File> <Optional: Source>'
  print 'INPUT GENEPRED FORMAT:'
  print 'genePred format that associates a gene name with gene prediction information'
  print 'No header on input (may support in the future)'
  print 'No comments on input (may support in future)'
  print '  ('
  print '  string  geneName;           "Name of gene as it appears in Genome Browser."'
  print '  string  name;               "Name of gene"'
  print '  string  chrom;              "Chromosome name"'
  print '  char[1] strand;             "+ or - for strand"'
  print '  uint    txStart;            "Transcription start position"'
  print '  uint    txEnd;              "Transcription end position"'
  print '  uint    cdsStart;           "Coding region start"'
  print '  uint    cdsEnd;             "Coding region end"'
  print '  uint    exonCount;          "Number of exons"'
  print '  uint[exonCount] exonStarts; "Exon start positions"'
  print '  uint[exonCount] exonEnds;   "Exon end positions"'
  print '  )'
  print 'OUTPUT GTF FORMAT:'
  print 'genePred format that associates a gene name with gene prediction information'
  print 'No header on output (may support in the future)'
  print 'No comments on output (may support in future)'
  print 'Frame is not implemented on output (may support in future)'
  print 'Source is an optional 2nd input'
  print '  ('
  print '  string  chrom;           "Chromsome name."'
  print '  string  source;          "Source e.g. hg19"'
  print '  string  feature;         "exon"'
  print '  uint    start;           "Feature start position"'
  print '  uint    end;             "Feature end position"'
  print '  float   score;           ". (empty)"'
  print '  char[1] strand;          "+ or - for strand"'
  print '  char[1] frame;           ". (empty) for exons"'
  print '  string  attributes;      "gene/transcript/exon ids"'
  print '  )'
  sys.exit()

linenum = 0
source = '.'
if(len(sys.argv) > 2):
  source = sys.argv[2]
with open(sys.argv[1]) as fp:
  for line in fp:
    linenum+=1
    if not line.startswith("#"):
      #if(len(line.split("\t")) != 11):
      #  print 'unexpected line length of '+str(len(line.split("\t")))+' on line '+str(linenum)
      #  sys.exit()
      vals  = line.rstrip("\r\n").split("\t")
      gene = vals[0]
      txn = vals[1]
      chrom = vals[2]
      strand = vals[3]
      txStart = vals[4]
      txEnd = vals[5]
      cdsStart = vals[6]
      cdsEnd = vals[7]
      exonCount = vals[8]
      exonStarts = vals[9]
      exonEnds = vals[10]
      exonStarts = exonStarts.rstrip(",").split(",")
      exonEnds = exonEnds.rstrip(",").split(",")
      for i in range (0,len(exonStarts)):
        exonnumber = i+1
        exonid = txn+'.'+str(exonnumber)
        print chrom + "\t" + source + "\texon\t" + str(int(exonStarts[i])+1) +"\t" + exonEnds[i] + "\t.\t" + strand + "\t.\t" + 'gene_id "'+gene+'"; transcript_id "'+txn+'"; exon_number "'+str(exonnumber)+'"; exon_id "'+exonid+'"; gene_name "'+gene+'";'
