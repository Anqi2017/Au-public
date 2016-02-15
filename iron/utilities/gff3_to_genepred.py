#!/usr/bin/python
import argparse, sys, re

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('infile',help="FILENAME if '-' then use STDIN")
  parser.add_argument('--transcript_only',action='store_true',help="use transcript name as genename")
  args = parser.parse_args()

  inf = sys.stdin
  if args.infile != '-':
    inf = open(args.infile)
  mRNA = {}
  exon = {}
  gene = {}
  # read in everything from gff3
  for line in inf:
    if re.match('^#',line): continue # skip to next if we're on a comment
    f = line.rstrip().split("\t")
    e = {}
    e['chrom'] = f[0]
    e['feature'] = f[2]
    e['start'] = int(f[3])  #1-base
    e['finish'] = int(f[4]) #1-base
    e['strand'] = f[6]
    group = f[8].rstrip(";").split(";")
    #print group
    for gone in group:
      if not re.search("=",gone): continue
      glist = gone.split("=")
      e[glist[0]] = glist[1]
    if e['feature'] == 'gene':
      if e['ID'] in gene:
        sys.stderr.write("WARNING duplicate entry of gene "+ e['ID']+"\n")
      gene[e['ID']] = e
    if e['feature'] == 'mRNA':
      if e['ID'] in mRNA:
        sys.stderr.write("WARNING duplicate entry of mRNA "+ e['ID']+"\n")
      mRNA[e['ID']] = e
    if e['feature'] == 'exon':
      if e['Parent'] not in exon:
        exon[e['Parent']] = {}
      if e['start'] in exon[e['Parent']]:
        sys.stderr.write("WARNING duplicate start of exon "+ json.dumps(e)+"\n")
      exon[e['Parent']][e['start']] = e
  # finished reading everything
  for tx in mRNA:  #for each transcript
    if not args.transcript_only: tx_parent = mRNA[tx]['Parent']
    if args.transcript_only: gid = tx
    else:
      gid = gene[tx_parent]['ID']
    starts = sorted(exon[tx].keys())
    sites = []
    chrom = False
    strand = False
    for start in starts:
      chrom = exon[tx][start]['chrom']
      strand = exon[tx][start]['strand']
      sites.append([exon[tx][start]['start']-1,exon[tx][start]['finish']])
    starts = [str(x[0]) for x in sites]
    ends = [str(x[1]) for x in sites]
    start = starts[0]
    finish = ends[len(ends)-1]
    print gid + "\t" + tx + "\t" + chrom + "\t" + strand + "\t" \
        + str(start) + "\t" + str(finish) + "\t" \
        + str(start) + "\t" + str(finish) + "\t" \
        + str(len(starts)) + "\t" + ",".join(starts)+",\t"+",".join(ends)+","

main()
