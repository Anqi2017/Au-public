#!/usr/bin/python
import sys, argparse
from GenePredBasics import GenePredDualLocusStream, GenePredEntry
import GenePredFuzzyBasics
from subprocess import Popen, PIPE
from multiprocessing import Lock, Pool, cpu_count

outind = 0
glock = Lock()

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('LR_sorted')
  parser.add_argument('SR_sorted')
  parser.add_argument('--fill_gaps',type=int,default=68)
  parser.add_argument('--required_evidence',type=int,default=1)
  parser.add_argument('--junction_tolerance',type=int,default=10)
  parser.add_argument('--threads',type=int,default=cpu_count())
  args = parser.parse_args()
  #for now lets assume LR psl and SR as bam

  # Get short stream ready
  cmd1 = "samtools view "+args.SR_sorted
  p1 = Popen(cmd1.split(),stdout=PIPE)
  cmd2 = "sam_to_psl.py -"
  p2 = Popen(cmd2.split(),stdin=p1.stdout,stdout=PIPE)
  cmd3 = "psl_to_target_genepred.py - --fill_gaps "+str(args.fill_gaps)
  p3 = Popen(cmd3.split(),stdin=p2.stdout,stdout=PIPE)
  shortstream = p3.stdout

  # Get long stream ready
  cmd4 = "psl_to_target_genepred.py "+args.LR_sorted+" --fill_gaps "+str(args.fill_gaps)
  p4 = Popen(cmd4.split(),stdout=PIPE)
  longstream = p4.stdout
  ds = GenePredDualLocusStream(longstream,shortstream)
  if args.threads > 1:
    p = Pool(processes=args.threads)
  while True:
    entry = ds.read_locus()
    if not entry: break
    if args.threads == 1:
      outdata = process_locus(entry,args)
      do_outs(outdata)
    else:
      p.apply_async(process_locus,args=(entry,args),callback=do_outs)
  if args.threads > 1:
    p.close()
    p.join()

def do_outs(outdata):
  if not outdata: return
  [outputs,totalrange] = outdata
  global glock
  glock.acquire()
  global outind
  sys.stderr.write(totalrange.get_range_string()+" "+str(outind)+"   \r")
  for output in outputs:
    print output
  glock.release()

def process_locus(entry,args):
  global glock
  global outind
  lr = entry[0]
  sr = entry[1]
  srjun = []
  for gpd in sr:
    if gpd.get_exon_count() > 1:
      srjun.append(gpd)
  if len(lr) == 0: return None
  totalrange = get_total_range(entry)
  #print '^^^^ Locus ^^^^'
  #print totalrange.get_range_string()
  #print str(len(lr))+"\t"+str(len(sr))+"\t"+str(len(srjun))
  # Get fuzzys from of all short reads
  srfzs = [GenePredFuzzyBasics.FuzzyGenePred(x) for x in srjun]
  for i in range(0,len(srfzs)): srfzs[i].gpds[0].entry['name'] = 'SR_'+str(i)
  fzs = GenePredFuzzyBasics.greedy_gpd_list_to_combined_fuzzy_list(lr,args.junction_tolerance)
  #print str(len(fzs)) + " genepreds"
  cnt = 0
  outind = 0
  outputs = []
  for fz in fzs:
    for i in range(0,len(fz.gpds)): 
      cnt += 1
      fz.gpds[0].entry['name'] = 'LR_'+str(cnt)
    g = GenePredEntry(fz.get_genepred_line())
    #print g.get_bed().get_range_string() + "\t" + str(g.get_exon_count())+" exons"
    parts = evaluate_junctions(fz,srfzs,args)
    for part in parts:
      glock.acquire()
      outind += 1
      glock.release()
      full = "LR_"+str(outind)+"\t"+"LR_"+str(outind)+"\t"+part
      outputs.append(full)
  return [outputs,totalrange]

def evaluate_junctions(fz,srfzs,args):
  cnt = 0
  working = fz.copy()
  if len(working.fuzzy_junctions) == 0: return []
  for i in range(0,len(working.fuzzy_junctions)):
    newjun = working.fuzzy_junctions[i]
    newjun.left.get_payload()['junc'] = []
    newjun.right.get_payload()['junc'] = []
    oldjun = fz.fuzzy_junctions[i]
    for srfz in srfzs:
      for sjun in srfz.fuzzy_junctions:
         if oldjun.overlaps(sjun,args.junction_tolerance):
           newjun.left.get_payload()['junc'].append(sjun.left.get_payload()['junc'][0])
           newjun.right.get_payload()['junc'].append(sjun.right.get_payload()['junc'][0])
           cnt +=1
  juncs = []
  starts = []
  ends = []
  evidences = []
  for i in range(0,len(fz.fuzzy_junctions)):
    evidence = len(working.fuzzy_junctions[i].left.get_payload()['junc'])
    if evidence >= args.required_evidence:
      if i == 0:
        starts.append(working.start.start+1)
      elif working.fuzzy_junctions[i].left.get_payload()['start']:
        starts.append(working.fuzzy_junctions[i].left.get_payload()['start'].start+1)
      else:
        starts.append(working.fuzzy_junctions[i-1].right.start+1)
      #now ends
      if i == len(fz.fuzzy_junctions)-1:
        ends.append(working.end.end)
      elif working.fuzzy_junctions[i].right.get_payload()['end']:
        ends.append(working.fuzzy_junctions[i].right.get_payload()['end'].end)
      else:
        ends.append(working.fuzzy_junctions[i+1].left.end)
      bestleft = GenePredFuzzyBasics.mode(working.fuzzy_junctions[i].left.get_payload()['junc'])
      bestright = GenePredFuzzyBasics.mode(working.fuzzy_junctions[i].right.get_payload()['junc'])
      juncs.append([bestleft,bestright])
      #print 'jun '+str(i)+' evid: '+str(evidence)+" "+str(bestleft)+" "+str(bestright)
    else:
      starts.append([])
      ends.append([])
      juncs.append([])
    evidences.append(evidence)
  #print juncs
  #print starts
  #print ends
  #print evidences
  # now we can put together the runs
  runs = []
  current_run = []
  for i in range(0,len(evidences)):
    if evidences[i] < args.required_evidence:
      if len(current_run) > 0:
        runs.append(current_run)
      current_run = []
      continue
    current_run.append(i)
  if len(current_run) > 0:
    runs.append(current_run)
  # now the runs are in runs
  #print 'runs:'
  parts = []
  for run in runs:
    sarr = []
    sarr.append(starts[run[0]]-1) #put back to zero index
    earr = []
    for i in range(0,len(run)):
      sarr.append(juncs[run[i]][1]-1)
      earr.append(juncs[run[i]][0])
    earr.append(ends[run[-1]])
    # ready to build a genepred!
    part = ''
    part += str(working.start.chr)+"\t"
    part += '+'+"\t"
    part += str(sarr[0])+"\t"
    part += str(earr[-1])+"\t"
    part += str(sarr[0])+"\t"
    part += str(earr[-1])+"\t"
    part += str(len(sarr))+"\t"
    part += ','.join([str(x) for x in sarr])+','+"\t"
    part += ','.join([str(x) for x in earr])+','
    parts.append(part)
  return parts

def get_total_range(entry):
  r = None
  if len(entry[0]) > 0:
    r = entry[0][0].get_bed()
    r.direction = None
  elif len(entry[1]) > 0:
    r = entry[1][0].get_bed()
    r.direction = None
  else:
    sys.stderr.write("ERROR\n")
    sys.exit()
  for e0 in entry[0]:
    b = e0.get_bed()
    b.direction = None
    r = r.merge(b)
  for e1 in entry[1]:
    b = e1.get_bed()
    b.direction = None
    r = r.merge(b)
  return r

if __name__=="__main__":
  main()
