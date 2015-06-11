#!/usr/bin/python
import SamBasics, SequenceBasics
from PSLBasics import line_to_entry as psl_line_to_entry
import argparse, sys, json

# This script should bridge BWA-mem sam output to another program that can 
# solve the splicing structure
# Pre: Take 1.  A BWA-mem sam file (unsorted raw output, otherwise it must be sorted by query name)
#           2.  The query reads from either --fasta or --fastq format
#           3.  The genomic reference fasta
# Post: Data where each row contains a json encoded data structure with a locus
#       to be solved by pairagon software.
#       It is an json encoded array containing [cdnaoffset, referenceoffset, seed fasta, cDNA fasta, reference fasta]
#       Where reference offset tells you how many bases your reference is offset
#

def main():
  parser = argparse.ArgumentParser(description="get guide from bwa-mem")
  parser.add_argument('sam_file',help="SAMFILE directly from bwa-mem or '-' for STDIN")
  parser.add_argument('--fasta',help="FASTAFILE query reads")
  parser.add_argument('--fastq',help="FASTQFILE query reads")
  parser.add_argument('--reference',help="FASTAFILE reference",required=True)
  args = parser.parse_args()
  inf = sys.stdin
  if args.sam_file != '-': inf = open(args.sam_file)
  reads = None
  reference = None
  if args.fasta:
    reads = SequenceBasics.read_fasta_into_hash(args.fasta)
  if args.reference:
    reference = SequenceBasics.read_fasta_into_hash(args.reference)
  if args.fastq:
    reads = {}
    gffr = SequenceBasics.GenericFastqFileReader(args.fastq)
    while True:
      e = gffr.read_entry()
      if not e: break
      reads[e['name']] = e['seq']
  if not reads:
    sys.stderr.write("Error reads are necessary\n")
    return
  spcf = SamBasics.SAMtoPSLconversionFactory()
  current_name = None
  locs = {}
  for line in inf:
    line = line.rstrip()
    if SamBasics.is_header(line):
      spcf.read_header_line(line)
      continue
    psl_line = spcf.convert_line(line)
    #print psl_line
    if not psl_line: continue
    e = psl_line_to_entry(psl_line)
    #print '----'
    #print "qstarts: " + str(e['qStarts'])
    #print "sizes: "+ str(e['blockSizes'])
    #print e['qStart']
    #print e['qEnd']
    name = e['qName']
    chrom = e['tName']
    for i in range(0,len(e['tStarts'])):
      loc = [e['tStarts'][i], e['tStarts'][i]+e['blockSizes'][i],e['qStarts'][i],e['qStarts'][i]+e['blockSizes'][i]]
      if name != current_name:
        # output buffer if its there
        if current_name:
          process_output(current_name,locs,reads,reference)
        locs = {}
        current_name = name
      if chrom not in locs:
        locs[chrom] = {}
      if e['strand'] not in locs[chrom]:
        locs[chrom][e['strand']] = []
      locs[chrom][e['strand']].append(loc)
    if current_name:
      process_output(current_name,locs,reads,reference)

def process_output(name,locs,reads,reference):
  #print name
  for chrom in locs:
    for strand in locs[chrom]:
      coords = locs[chrom][strand]
      # need to go through the coordinates and get loci
      scoords = sorted(coords,key=lambda x: x[0])
      #print chrom + " " + strand
      cset = []
      farthest = 0
      for c in scoords:
        if c[1] > farthest+1000000:
          if len(cset) > 0:  
            process_output2(name,chrom,strand,cset,reads,reference)
          cset = []
        cset.append(c)
        farthest = c[1]
      if len(cset) > 0:
        process_output2(name,chrom,strand,cset,reads,reference)

def process_output2(name,chrom,strand,coords,reads,reference):
  #coords = make_ungapped(coords)
  startpoint = coords[0][0]
  cdnastartpoint = coords[0][2]
  endpoint = coords[len(coords)-1][1]
  locus = reference[chrom][startpoint:endpoint]
  report = ''
  report += '>' + name + "\n"
  report += "genomic_boundary_start=1 genomic_boundary_end="+str(len(locus)) + " strand="+strand + "\n"
  report += "count="+str(len(coords)) + "\n"
  for c in coords:
    #print c
    report += '('+str(c[0]+1-startpoint)+', '+str(c[2]+1-cdnastartpoint)+') ('+str(c[1]-startpoint)+', '+str(c[3]-cdnastartpoint)+')' + "\n"
    #if reads:
    #  oread = reads[name]
    #  if strand == '-': oread = SequenceBasics.rc(oread)
    #  print oread[c[2]:c[3]].upper()
  #print report
  cDNAfasta = ">"+name+"\n"+reads[name][cdnastartpoint:]+"\n"
  targetfasta = ">genomic"+"\n"+locus.upper()+"\n"
  data = [cdnastartpoint, startpoint, report, cDNAfasta, targetfasta]
  print json.dumps(data)

def make_ungapped(coords):
  ungapped = []
  if len(coords) == 1: return coords
  ungapped.append(coords[0][:])
  for i in range(1,len(coords)):
    gapsize = abs(coords[i-1][3]-coords[i][2])
    ungapped.append([coords[i][0],coords[i][1],coords[i-1][3],coords[i][3]])
  return ungapped

if __name__=="__main__":
  main()
