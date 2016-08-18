#!/usr/bin/python
import sys, argparse, random
from GenePredBasics import GenePredDualLocusStream, GenePredEntry
import GenePredFuzzyBasics
from subprocess import Popen, PIPE
from multiprocessing import Lock, Pool, cpu_count

glock = Lock()
locus_count = 0
of = sys.stdout
of_table = None
def main():
  global of
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('LR_sorted',help="sorted gpd")
  parser.add_argument('SR_sorted',help="sorted short reads. be carefule these are alphabetical chromosome. samtools sort does not garuntee this.")
  parser.add_argument('-r','--reference',required=True,help="reference fasta")
  parser.add_argument('-o','--output',help="Leave unset for STDOUT")
  parser.add_argument('--fill_gaps',type=int,default=68,help="for psl to gpd")
  parser.add_argument('--required_evidence',type=int,default=1)
  parser.add_argument('-j','--junction_tolerance',type=int,default=10)
  parser.add_argument('--threads',type=int,default=cpu_count())
  #group = parser.add_mutually_exclusive_group()
  #group.add_argument('--by_locus',action='store_true')
  #group.add_argument('--by_read')
  #group1 = parser.add_mutually_exclusive_group()
  #group1.add_argument('--LR_bam')
  #group1.add_argument('--LR_psl',action='store_true')
  #group1.add_argument('--LR_gpd',action='store_true')
  parser.add_argument('--downsample',type=int,default=500,help="Set to -1 for all reads, otherwise this is max read depth") #maximum short reads to consider at a junction site
  parser.add_argument('--maxlocusSR',type=int,default=50000,help="Maximum number of short reads for locus") 
  parser.add_argument('--output_original_table',help="help trace results back to original")
  args = parser.parse_args()
  if args.output: of = open(args.output,'w')
  if args.output_original_table:
    global of_table
    of_table = open(args.output_original_table,'w')

  #for now lets assume LR psl and SR as bam
  # Get short stream ready
  cmd1 = "samtools view -F 4 "+args.SR_sorted
  sys.stderr.write(cmd1+"\n")
  p1 = Popen(cmd1.split(),stdout=PIPE)
  cmd2 = "awk '$6~/N/'"
  p2 = Popen(cmd2,stdin=p1.stdout,stdout=PIPE,shell=True)
  # we never need more than the required evidence at a particular position
  # more non positional duplicates can still help choose the best when many
  # reads are available
  cmd3 = "filter_sam_positional_duplicates.py --positional_duplicates "+str(args.required_evidence)+" - "
  p3 = Popen(cmd3.split(),stdin=p2.stdout,stdout=PIPE)
  cmd4 = "sam_to_psl.py - -r "+args.reference
  p4 = Popen(cmd4.split(),stdin=p3.stdout,stdout=PIPE)
  cmd5 = "psl_to_target_genepred.py - --fill_gaps "+str(args.fill_gaps)
  p5 = Popen(cmd5.split(),stdin=p4.stdout,stdout=PIPE)
  shortstream = p5.stdout

  # Get long stream ready
  longstream = None
  #if not args.LR_gpd:
  #  cmd5 = "psl_to_target_genepred.py "+args.LR_sorted+" --fill_gaps "+str(args.fill_gaps)
  #  p5 = Popen(cmd5.split(),stdout=PIPE)
  #  longstream = p5.stdout
  #else:  # we have gpd
  longstream = open(args.LR_sorted)
  ds = GenePredDualLocusStream(longstream,shortstream)
  #p = None
  if args.threads > 1:
    p = Pool(processes=args.threads)
  while True:
    entry = ds.read_locus()
    if not entry: break
    # lets clean up the short reads here a bit
    if len(entry[0]) == 0: continue
    sr = []
    for srgpd in entry[1]:
      if srgpd.get_exon_count() < 2: continue
      sr.append(srgpd)
    if len(sr) == 0: continue
    if len(sr) > args.maxlocusSR:
      sys.stderr.write("\nWARNING: max locus SR exceeded "+sr[0].get_bed().get_range_string()+" with "+str(len(sr))+"\n")
      newsr = sr
      random.shuffle(newsr)
      sr = newsr[0:args.maxlocusSR]
      sys.stderr.write("\n reduced to "+str(len(sr))+"\n")
    if args.threads == 1:
      outdata = process_locus(entry[0],sr,args)
      do_outs(outdata)
      #sys.stderr.write("\nfinished locus\n")
    else:
      p.apply_async(process_locus,args=(entry[0],sr,args),callback=do_outs)
  if args.threads > 1:
    p.close()
    p.join()

def do_outs(outdata):
  if not outdata: return
  if len(outdata)==0: return
  outputs,totalrange = outdata
  global locus_count
  global glock
  global of
  global of_table
  glock.acquire()
  for output_combo in outputs:
    output, source_names = output_combo
    locus_count += 1
    of.write('LR_'+str(locus_count)+"\t"+'LR_'+str(locus_count)+"\t"+output+"\n")
    if of_table:
      for name in source_names:
        of_table.write('LR_'+str(locus_count)+"\t"+name+"\n")
        
  sys.stderr.write(totalrange.get_range_string()+" "+str(locus_count)+"        \r")
  glock.release()

def process_locus(lr,srin,args):
  if len(lr) == 0: return None
  totalrange = get_total_range(lr)
  #print '^^^^ Locus ^^^^'
  #print totalrange.get_range_string()
  #print str(len(lr))+"\t"+str(len(sr))+"\t"+str(len(srjun))
  # Get fuzzys from of all short reads
  sr = {}
  #do this more time consuming cutdown ont he SR data after sending to a thread
  for srgpd in srin:
      srfz = GenePredFuzzyBasics.FuzzyGenePred(srgpd,juntol=args.junction_tolerance)
      for j in srfz.fuzzy_junctions:
        junstr = j.left.chr+':'+str(j.left.end)+','+str(j.right.end)
        if junstr not in sr:
          sr[junstr] = {}
          sr[junstr]['cnt'] = 0
          sr[junstr]['fzjun'] = j
        sr[junstr]['cnt'] += 1

  #srfzs = [GenePredFuzzyBasics.FuzzyGenePred(x) for x in srjun]
  #for i in range(0,len(srfzs)): srfzs[i].gpds[0].entry['name'] = 'SR_'+str(i)
  fzs = GenePredFuzzyBasics.greedy_gpd_list_to_combined_fuzzy_list(lr,args.junction_tolerance)
  #print str(len(fzs)) + " genepreds"
  outputs = []
  #if args.threads > 1:
  #  p = Pool(processes=args.threads)
  for fz in fzs:
    #if args.by_read:
    #  if args.threads > 1 and args.by_read:
    #    p.apply_async(do_fuzzy,args=(fz,sr,args),callback=do_outs)
    #  else:
    #    outs = do_fuzzy(fz,sr,args)
    #    do_outs([outs,totalrange])
    #else:
    outs = do_fuzzy(fz,sr,args)
    #  do_outs([outs,totalrange])
    for o in outs: outputs.append(o)
  #if args.threads > 1 and args.by_read:
  #  p.close()
  #  p.join()
  #if not args.by_read:
  return [outputs,totalrange]
  #return

def do_fuzzy(fz,sr,args):
    outputs = []
    cnt = 0
    for i in range(0,len(fz.gpds)): 
      cnt += 1
      #fz.gpds[0].entry['name'] = 'LR_'+str(cnt)
    g = GenePredEntry(fz.get_genepred_line())
    #print g.get_bed().get_range_string() + "\t" + str(g.get_exon_count())+" exons"
    parts = evaluate_junctions(fz,sr,args)
    for part in parts:
      #full = "LR_"+str(outind)+"\t"+"LR_"+str(outind)+"\t"+part
      outputs.append(part)
    return outputs

def evaluate_junctions(fz,sr,args):
  cnt = 0
  source_names = [x.entry['name'] for x in fz.gpds]
  working = fz.copy()
  if len(working.fuzzy_junctions) == 0: return []
  for i in range(0,len(working.fuzzy_junctions)):
    newjun = working.fuzzy_junctions[i]
    newjun.left.get_payload()['junc'] = []
    newjun.right.get_payload()['junc'] = []
    oldjun = fz.fuzzy_junctions[i]
    for srjun in sr:
         sjun = sr[srjun]['fzjun']
         if oldjun.overlaps(sjun,args.junction_tolerance):
           for i in range(0,min(sr[srjun]['cnt'],args.downsample)):
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
        starts.append(working.start.start)
      elif working.fuzzy_junctions[i].left.get_payload()['start']:
        starts.append(working.fuzzy_junctions[i].left.get_payload()['start'].start)
      else:
        starts.append(working.fuzzy_junctions[i-1].right.start)
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
    # Final quality check here
    gpd = GenePredEntry("test1\ttest1\t"+part)
    if not gpd.is_valid():
      sys.stderr.write("\nWARNING skipping invalid GPD\n"+gpd.get_line()+"\n")
      continue
    parts.append([part,source_names])
  #print parts
  return parts

def get_total_range(entry):
  r = None
  if len(entry) > 0:
    r = entry[0].get_bed()
    r.direction = None
  else:
    sys.stderr.write("ERROR\n")
    sys.exit()
  for e0 in entry:
    b = e0.get_bed()
    b.direction = None
    r = r.merge(b)
  return r

if __name__=="__main__":
  main()
