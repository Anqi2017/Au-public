#!/usr/bin/python
import sys, argparse, itertools
from Bio.Format.Sam import BAMFile
from Bio.Range import GenomicRange

def main():
  parser = argparse.ArgumentParser(description="Take a BAM file with multiple paths and determine gapped alignments",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="BAM file")
  parser.add_argument('--max_query_overlap',type=int,default=10,help="Consider two alignments incompatible if greater than this")
  parser.add_argument('--max_target_overlap',type=int,default=10,help="Consider two alignments incompatible if greater than this")
  parser.add_argument('--max_target_gap',type=int,default=500000,help="Consider a gapped alignment incompatible if greater than this")
  parser.add_argument('--max_query_gap',type=int,help="Consider a gapped alignment incompatible if greater thant this")
  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="Result should be this much better than the original")
  parser.add_argument('--output','-o',help="Output file or stdout if not defined")
  args = parser.parse_args()

  #query_lengths = {}
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  bf = BAMFile(args.input,skip_index=True)
  #use our memory
  aligned = {}
  unaligned = {}
  z = 0
  for e in bf:
    z+=1
    if z%1000==0: sys.stderr.write(str(z)+"\r")
    if e.is_aligned():
      #q = e.get_query_range()
      o = e.get_original_query_length()
      oarng = e.get_actual_original_query_range()
      b = e.get_aligned_bases_count()
      t = e.get_target_range()
      v = {'dir':oarng.direction,'qlen':o,'oarng':oarng,'trng':t,'abases':b}
      if oarng.chr not in aligned: aligned[oarng.chr] = []
      aligned[oarng.chr].append(v)
    else: 
      unaligned[e.value('qname')] = e.get_query_length()
  sys.stderr.write("\n")
  for qname in aligned:
    paths = get_compatible_paths(aligned[qname],args)
    #print len(paths)
    #print sum([x['abases'] for x in paths])
    rngs = [x['oarng'] for x in sorted(paths['overall'],key=lambda x: x['oarng'].start)]
    tot = [rngs[0]]
    if len(rngs) > 1:
      for r in rngs[1:]:
        if r.overlaps(tot[-1]): tot[-1] = tot[-1].merge(r)
        else: tot.append(r)
    bcount = sum([x.length() for x in tot])
    scount = paths['single']['oarng'].length()
    #print bcount
    of.write(str(len(paths['overall']))+"\t"+str(scount)+"\t"+str(bcount)+"\t"+str(paths['overall'][0]['qlen'])+"\n")
  for qname in unaligned:
    of.write("0\t0\t0\t"+str(unaligned[qname])+"\n")

# given the alignments return the compatible paths
def get_compatible_paths(alns,args):
  possible = {}
  for a in alns:
    if a['trng'].chr not in possible: possible[a['trng'].chr] = {}
    if a['dir'] not in possible[a['trng'].chr]: possible[a['trng'].chr][a['dir']] = []
    possible[a['trng'].chr][a['dir']].append(a)
  bestlone = None
  bestoutsingle = None
  for chr in possible:
    for dir in possible[chr]:
      poss = possible[chr][dir]
      (lone,curr) = sorted([[x['abases'],x] for x in poss])[-1]
      if not bestlone:
        bestlone = lone
        bestoutsingle = curr
      if lone > bestlone:
        bestlone = lone
        bestoutsingle = curr
  output = None
  bestoutput  = None
  bestbases = 0
  for chr in possible:
    for dir in possible[chr]:
      poss = possible[chr][dir]
      #print str(len(poss))+' members'
      isets = get_index_sets(len(poss))
      for subset in isets:
        v = analyze_subset([poss[i] for i in subset],args)
        if v:
          effective = v['aligned'] - v['overlap']
          if float(effective)>float(bestlone)+ float(bestlone)*args.required_fractional_improvement and effective > bestbases:
            bestoutput = [poss[i] for i in subset]
            bestbases = effective
          #print effective
  if bestoutput: return {'single':bestoutsingle,'overall':bestoutput}
  return {'single':bestoutsingle,'overall':[bestoutsingle]}

def analyze_subset(subset,args):
  t_overlap = 0
  q_overlap = 0
  total_q_overlap = 0
  total_abases = sum([x['abases'] for x in subset])
  for i in range(len(subset)):
    for j in range(i+1,len(subset)):
      ov = subset[i]['oarng'].overlap_size(subset[j]['oarng'])
      if ov > q_overlap: q_overlap = ov
      tv = subset[i]['trng'].overlap_size(subset[j]['trng'])
      if tv > t_overlap: t_overlap = ov
      total_q_overlap += ov
  #print t_overlap
  #print q_overlap
  if t_overlap > args.max_target_overlap: return False
  if q_overlap > args.max_query_overlap: return False
  tordered = sorted(subset,key=lambda x: x['trng'].start)
  for i in range(len(tordered)-1):
    v = tordered[i]
    vnext = tordered[i+1]
    t = v['trng']
    q = v['oarng']
    tnext = vnext['trng']
    qnext = vnext['oarng']
    tgap = tnext.start - t.end-1
    if tgap > args.max_target_gap:
      return False
    #print t.get_range_string()+"    "+tnext.get_range_string()
    #print q.get_range_string()+"    "+qnext.get_range_string()
    qgap = qnext.start - q.end - 1
    if v['dir'] == '-': qgap = q.start-qnext.end-1
    #print qgap
    if qgap < args.max_query_overlap*-1: return False
    if args.max_query_gap:
      if qgap > args.max_query_gap: return False
    ###### if we are still here, this is okay #######
  return {'overlap':total_q_overlap,'aligned':total_abases}

def get_index_sets(indlen):
  r = []
  inds = range(0,indlen)
  for l in range(1,len(inds)+1):
    for subset in itertools.combinations(inds,l):
      r.append(subset)
  return r

if __name__=="__main__":
  main()
