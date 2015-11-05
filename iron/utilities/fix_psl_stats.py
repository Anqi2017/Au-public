#!/usr/bin/python
import sys, argparse
import PSLBasics
from SequenceBasics import read_fasta_into_hash, FastaHandleReader, rc
# Fixes the matches/mismatches/Ncount of a PSL file
# This file will take as an input, a PSL file,
# a query file, and a reference file.
# Pre: Every query sequence is required to be present
#      in the query file and every reference in the reference file
# IMPORTANT!!!! input PSL is ordered by query, and query FASTA is ordered by query

def main():
  parser = argparse.ArgumentParser(description="Correct the matches/mismatches and Ncount of a PSL file")
  parser.add_argument('input',help="PSLFILE or - for STIDN")
  parser.add_argument('reference',help="FASTAFILE reference genome")
  parser.add_argument('query',help="FASTAFILE query sequences")
  parser.add_argument('--minimum_intron_size',type=int,default=68,help="INT")
  parser.add_argument('--ordered_query',action='store_true',help="The query psl and fasta are both ordered by query name for optimal performance")
  args = parser.parse_args()
  # Read in the reference genome
  sys.stderr.write("Reading in reference genome\n")
  g = read_fasta_into_hash(args.reference)
  sys.stderr.write("Finished reading "+str(len(g.keys()))+" reference sequences\n")
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  fhr = FastaHandleReader(open(args.query))
  last_fasta = fhr.read_entry()
  if not last_fasta:
    sys.stderr.write("ERROR: No query sequences\n")
    sys.exit()
  for line in inf:
    p = PSLBasics.PSL(line)
    if not p.validate():
      sys.stderr.write("WARNING skipping invalid PSL entry. This script fixes improper mismatch and match counts .. and gap counts... it doesn't perform miracles. Problem line:\n"+line.rstrip()+"\n")
    n = p.value('qName')
    if not last_fasta:
      sys.stderr.write("ERROR: Ran out of query sequences too soon.  Are they sorted properly\n")
      sys.exit()
    while last_fasta['name'] != n:
      last_fasta = fhr.read_entry()
    f = last_fasta
    nCount = 0
    matches = 0
    misMatches = 0
    prev_qE = 0
    prev_tE = 0
    qNumInsert = 0
    qBaseInsert = 0
    tNumInsert = 0
    tBaseInsert = 0
    for i in range(p.value('blockCount')):
      blen = p.value('blockSizes')[i]
      qS = p.value('qStarts')[i] #query start
      qE = qS + blen             #query end
      tS = p.value('tStarts')[i] #target start
      tE = tS + blen             #target end
      #Work on gaps
      if prev_qE > 0 or prev_tE > 0: #if its not our first time through
        tgap = tS-prev_tE
        if tgap < args.minimum_intron_size and tgap > 0:
          tNumInsert += 1
          tBaseInsert += tgap
        qgap = qS-prev_qE
        if qgap > 0:
          qNumInsert += 1
          qBaseInsert += qgap
      query = f['seq']
      if p.value('strand') == '-':
        query = rc(f['seq'])
      qseq = query[qS:qE].upper()
      rseq = g[p.value('tName')][tS:tE].upper()
      #print qseq+"\n"+rseq+"\n"
      for j in range(0,blen):
        if qseq[j] == 'N':
          nCount += 1
        elif qseq[j] == rseq[j]:
          matches += 1
        else:
          misMatches += 1
      prev_qE = qE
      prev_tE = tE
    p.entry['matches'] = matches
    p.entry['misMatches'] = misMatches
    p.entry['nCount'] = nCount
    p.entry['qNumInsert'] = qNumInsert
    p.entry['qBaseInsert'] = qBaseInsert
    p.entry['tNumInsert'] = tNumInsert
    p.entry['tBaseInsert'] = tBaseInsert
    p.entry['qSize'] = len(query)
    p.entry['tSize'] = len(g[p.value('tName')]) 
    print p.get_line()
  fhr.close()
if __name__=="__main__":
  main()
