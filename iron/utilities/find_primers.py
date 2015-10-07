#!/usr/bin/python
import argparse, sys, re
from SequenceBasics import read_fasta_into_hash
from PairwiseAlignmentBasics import SmithWatermanAligner

def main():
  parser = argparse.ArgumentParser(description="Find frequency of poly-A tails in genomic sequences")
  parser.add_argument('transcriptome_fasta',help="FASTA_FILE genome")
  args = parser.parse_args()
  tx = read_fasta_into_hash(args.transcriptome_fasta)
  prog1 = re.compile('(A+)')
  prog2 = re.compile('(T+)')
  totals = {}
  total_length = 0
  lenlow = 18
  lenhigh = 18
  counts = {}
  z = 0
  for name in tx:
    z+=1
    sys.stderr.write(str(z)+"\r")
    seq = tx[name].upper()
    explode(counts,seq,lenlow,lenhigh)  
    longest = 0
    if z %20 == 0: 
      for part in counts.keys():
        if counts[part] <= 3: del counts[part]
  sys.stderr.write("\n")
  numuse = 200
  rankedsets= {}
  z = 0
  for myset in [[x,counts[x]] for x in sorted(counts, key=counts.get,reverse=True)][0:numuse]:
    rankedsets[z] = myset
    z +=1
  lastsets = -1
  numsets = len(rankedsets.keys())
  while numsets != lastsets:
    sys.stderr.write(str(numsets)+"          \r")
    lastsets = numsets
    reduceranks(rankedsets) # remove highly similar sets
    numsets = len(rankedsets.keys())
  sys.stderr.write("\n")
  for i in sorted(rankedsets.keys()):
    print rankedsets[i][0]+"\t"+str(rankedsets[i][1])

def reduceranks(rankedsets):
  rankspresent = sorted(rankedsets.keys())
  #print rankspresent
  for i in range(0,len(rankspresent)):
    irank = rankspresent[i]
    for j in range(i+1,len(rankspresent)):
      jrank = rankspresent[j]
      swa = SmithWatermanAligner()
      swa.set_sequences(rankedsets[irank][0],rankedsets[jrank][0])
      alignment = swa.align()
      target_length = len(rankedsets[irank][0])
      matches = alignment.count_matches()
      if matches >= target_length-2 and (alignment.count_misMatches()+alignment.count_qBaseInsert()+alignment.count_tBaseInsert()) <= 2:
        #print rankedsets[irank]
        #print rankedsets[jrank]
        #alignment.print_alignment()
        rankedsets[irank][1] += rankedsets[jrank][1]
        del rankedsets[jrank]
        #print '---'
        return

def explode(counts,seq,lenlow,lenhigh):
  for i in range(lenlow,lenhigh+1):
    for j in range(0,len(seq)-i):
      part = seq[j:j+i]
      if part not in counts: counts[part] = 0
      counts[part] += 1
  return

if __name__=="__main__":
  main()
