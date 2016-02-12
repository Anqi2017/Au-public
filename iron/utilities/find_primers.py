#!/usr/bin/python
import argparse, sys, re
from SequenceBasics import FastaHandleReader, rc
from PairwiseAlignmentBasics import SmithWatermanAligner

def main():
  parser = argparse.ArgumentParser(description="Find primers in a sequence")
  parser.add_argument('input',help="FASTA_FILE genome or - for STDIN")
  parser.add_argument('--AT_end_limit',type=int,default=4,help='Maxmimum number of A/T to look for at the end')
  parser.add_argument('--overlap_join',type=int,default=8,help='Join together matches with this much exact overlap')
  parser.add_argument('--end_criteria',type=int,default=5000,help='Stop when you have seen a k-mer this many times')
  parser.add_argument('--total_candidates',type=int,default=100,help='Look at this number of candidates')
  parser.add_argument('--kmersize',type=int,default=18,help='Look at this number of candidates')
  args = parser.parse_args()

  if args.input == '-':
    args.input = sys.stdin
  else:
    args.input = open(args.input)

  #tx = read_fasta_into_hash(args.transcriptome_fasta)
  totals = {}
  total_length = 0
  lenlow = args.kmersize
  lenhigh = args.kmersize
  counts = {}
  z = 0
  reader = FastaHandleReader(args.input)
  while True:
    e = reader.read_entry()
    if not e: break
    z+=1
    sys.stderr.write(str(z)+"\r")
    seq = e['seq']
    explode(counts,seq,lenlow,lenhigh)  
    longest = 0
    if z %20 == 0: 
      for part in counts.keys():
        if counts[part] <= 3: del counts[part]
        else:
          edgemax = edgeAT(part)
          if edgemax > args.AT_end_limit: del counts[part]
      biggest =  sorted(counts, key=counts.get,reverse=True)
      if len(biggest) > 1:
        if counts[biggest[1]] > args.end_criteria:
          break
  sys.stderr.write("\n")
  numuse = args.total_candidates
  rankedsets= {}
  z = 0
  #for seq in counts.keys():
  #  edgemax = edgeAT(seq)
  #  if edgemax > args.AT_end_limit: del counts[seq]
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
  # now we have our best candidates
  #for i in sorted(rankedsets.keys()):
  #  print rankedsets[i][0]+"\t"+str(rankedsets[i][1])

  lastsets = 0
  numsets = len(rankedsets.keys())
  while numsets != lastsets:
    lastsets = numsets
    combine_overlapping(rankedsets,args)
    numsets = len(rankedsets.keys())
  for i in sorted(rankedsets.keys()):
    print rankedsets[i][0]+"\t"+str(rankedsets[i][1])

def combine_overlapping(rankedsets,args):
  rnums = rankedsets.keys()
  for i in range(0,len(rnums)):
    for j in range(i+1,len(rnums)):
      ov = overlapped(rankedsets[rnums[i]][0],rankedsets[rnums[j]][0],args.overlap_join)
      if ov:
        rankedsets[rnums[i]][1] += rankedsets[rnums[j]][1]
        rankedsets[rnums[i]][0] = ov
        del rankedsets[rnums[j]]
        return

def overlapped(seq1,seq2,thresh):
  v1 = starts_with(seq1,seq2,thresh)
  if v1:  return v1
  v2 = starts_with(seq2,seq1,thresh)
  if v2: return v2
  v3 = starts_with(seq1,rc(seq2),thresh)
  if v3: v3
  v4 = starts_with(seq2,rc(seq1),thresh)
  if v4: v4
  return False
   
# case of overlap where a sequence begins with seq1 then seq2 is continuing it
def starts_with(seq1,seq2,thresh):
  p1 = re.compile(seq2[0:thresh])
  m = p1.search(seq1)
  if m:
    #print 'over'
    #print seq1+"\t"+seq2
    #print m.start()
    total_fragment = seq1[m.start():]
    if re.match(total_fragment,seq2):
      # really have a match
      return seq1[0:m.start()]+seq2
  return False

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
        #use the set with the lower edgeAT value
        name1 = rankedsets[irank][0]
        name2 = rankedsets[jrank][0]
        if edgeAT(name1) > edgeAT(name2): rankedsets[irank][0] = name2 #switch to name 2
        #print '---'
        del rankedsets[jrank]
        return

def explode(counts,seq,lenlow,lenhigh):
  for i in range(lenlow,lenhigh+1):
    for j in range(0,len(seq)-i):
      part = seq[j:j+i]
      if part not in counts: counts[part] = 0
      counts[part] += 1
  return

def edgeAT(seq):
    prog1 = re.compile('^(A+)')
    prog2 = re.compile('^(T+)')  
    prog3 = re.compile('(A+)$')
    prog4 = re.compile('(T+)$')
    m1 = prog1.search(seq)
    m2 = prog2.search(seq)
    m3 = prog3.search(seq)
    m4 = prog4.search(seq)
    if m1: m1 = len(m1.group(1))
    else: m1 = 0
    if m2: m2 = len(m2.group(1))
    else: m2 = 0
    if m3: m3 = len(m3.group(1))
    else: m3 = 0
    if m4: m4 = len(m4.group(1))
    else: m4 = 0
    edgemax = max([m1,m2,m3,m4])
    return edgemax
if __name__=="__main__":
  main()
