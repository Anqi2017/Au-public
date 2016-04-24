#!/usr/bin/python
import argparse, sys
from RangeBasics import Loci, Locus
from MultiplePSLBasics import GenericOrderedMultipleAlignmentPSLReader as GOMAPR
from multiprocessing import Lock, Pool, cpu_count

best_report_fh = None
output_fh = None
glock = Lock()

# PRE: A PSL file that is sorted by query name
# POST: A report on gapped alignments and segments of the alignment contributing to a most likely path
#       A psl file with best paths listed by query name ordered by query start (not reverse complemented)

def main():
  global best_report_fh
  global output_fh
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="INPUT name sorted PSL")
  parser.add_argument('--minimum_alignment_coverage',type=int,default=1)
  parser.add_argument('--maximum_paths',type=int,help="maximum number of paths to consider")
  parser.add_argument('--multipath_score_improvement',type=float,default=0,help="Require this fraction of a score imporvement over the single score")
  group1 = parser.add_mutually_exclusive_group()
  group1.add_argument('--fusion',action='store_true')
  group1.add_argument('--maximum_intron',type=int,default=400000,help="the maximum distance to try to connect gapped alignments")
  parser.add_argument('--maximum_query_gap',type=int,default=-1)
  parser.add_argument('--maximum_query_overlap',type=int,default=20)
  parser.add_argument('--maximum_target_overlap',type=int,default=0)
  parser.add_argument('--maximum_query_fraction_overlap',type=float,default=0.2)
  parser.add_argument('--best_report',help="output file to save a report on the paths")
  parser.add_argument('-o','--output',default='-',help="default STDOUT output file for psl of best paths where each name has only entries contributing to the best path ordered by their query coordiantes.")
  parser.add_argument('--threads',default=cpu_count(),type=int,help="number of threads defautl cpu_count")
  args = parser.parse_args()
  if args.input=='-': args.input = sys.stdin
  else: args.input = open(args.input)
  if args.output=='-': output_fh = sys.stdout
  else: output_fh = open(args.output,'w')
  if args.best_report:
    best_report_fh = open(args.best_report,'w')
  seen_reads = set()
  mpslr = GOMAPR(args.input)
  if args.threads > 1:
    p = Pool(processes=args.threads)
  rcnt = 0
  while True:
    mpa = mpslr.read_next()
    if not mpa: break
    read_name = mpa.entries[0].value('qName')
    if read_name in seen_reads:
      sys.stderr.write("ERROR: input reads need to be ordered by read name\n")
      sys.exit()
    rcnt+=1
    if rcnt%1000==0: sys.stderr.write(str(rcnt)+"        \r")
    seen_reads.add(read_name)
    # Enforce minimum alignment coverage
    if args.minimum_alignment_coverage:
      newvals = []
      for e in mpa.entries:
        if e.get_coverage() >= args.minimum_alignment_coverage:
          newvals.append(e)
      mpa.entries = newvals
    # Restrict ourselves to the paths with the most coverage
    if args.maximum_paths:
      if len(mpa.entries) > args.maximum_paths:
        bestentries = sorted(mpa.entries[:], key=lambda x: x.get_coverage(),reverse=True)
        mpa.entries = bestentries[0:args.maximum_paths]
    if args.threads < 2:
      v = process_read(mpa,args)
      do_outputs(v)
    else:
      p.apply_async(process_read,args=(mpa,args),callback=do_outputs)
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\n")

def do_outputs(invals):
  if not invals: return
  report, psl_lines = invals
  global best_report_fh
  global output_fh
  global glock
  glock.acquire()
  if best_report_fh:
    best_report_fh.write(report+"\n")
  for line in psl_lines:
    output_fh.write(line+"\n")
  glock.release()  

def process_read(mpa,args):
    # Filter entries by a minimum alignment coverage
    newentries = []
    for i in [i for i in range(0,len(mpa.entries)) if mpa.entries[i].get_coverage() > args.minimum_alignment_coverage]:
      newentries.append(mpa.entries[i])
    mpa.entries = newentries

    # Find best singles
    bestsingle = None
    bestsinglescore = -1
    for i in range(0,len(mpa.entries)):
      totalcov = mpa.entries[i].get_coverage()
      weightedcov = float(mpa.entries[i].get_coverage())*float(mpa.entries[i].get_quality())
      if weightedcov > bestsinglescore:
        bestsinglescore = weightedcov
        bestsingle = i
    if bestsinglescore == -1: 
      sys.stderr.write("failed to find a single path\n")
      return None
    my_max_intron = args.maximum_intron
    if args.fusion: my_max_intron = -1 # we can look any distance for a group
    mpa.compatible_graph(max_intron=my_max_intron,max_query_overlap=args.maximum_query_overlap,max_gap=args.maximum_query_gap,max_target_overlap=args.maximum_target_overlap,max_query_fraction_overlap=args.maximum_query_fraction_overlap)
    ps = mpa.get_root_paths()
    bestpath = [bestsingle]
    bestscore = 0
    besttotalcov = 0
    allscores = []
    allcov = []
    best_path_index = -1
    zz = 0
    for path in ps:
      totalcov = sum([mpa.entries[i].get_coverage() for i in path])
      weightedcov = sum([float(mpa.entries[i].get_coverage())*float(mpa.entries[i].get_quality()) for i in path])
      allscores.append(weightedcov)
      allcov.append(totalcov)
      if weightedcov > bestscore: 
        bestscore = weightedcov
        bestpath = path
        besttotalcov = totalcov
        best_path_index = zz
      zz+=1
    #if not bestpath: return None
    otherpaths = []
    for i in range(0,len(ps)):
      if i != best_path_index:
        otherpaths.append(ps[i])
    query_target_coverages = []
    for other_path in otherpaths:
      qcov = 0
      tcov = 0
      for other_entry in [mpa.entries[i] for i in other_path]:
        for entry in [mpa.entries[j] for j in bestpath]:
          qcov += other_entry.query_overlap_size(entry)
          tcov += other_entry.target_overlap_size(entry)
      query_target_coverages.append(str(qcov)+'/'+str(tcov))

    gapsizes = []
    if len(bestpath) > 1:
      gapsizes = [mpa.entries[bestpath[j+1]].get_query_bed().start - mpa.entries[bestpath[j]].get_query_bed().end -1 for j in range(0,len(bestpath)-1)]
    #print mpa.g.get_status_string()
    #print [mpa.entries[i].get_target_bed().get_range_string() for i in bestpath]
    #print [mpa.entries[i].get_query_bed().get_range_string() for i in bestpath]
    #print [mpa.entries[i].get_quality() for i in bestpath]
    #print [mpa.entries[i].get_coverage() for i in bestpath]
    #print gapsizes
    #print bestscore
    #print bestsinglescore

    #See if we should use the single path score instead
    if len(path) > 1 and bestsinglescore*(1+args.multipath_score_improvement) > bestscore:
      bestpath = [bestsingle]
      besttotalcov = mpa.entries[bestsingle].get_coverage()
      bestscore = bestsinglescore
    query_span = mpa.entries[bestpath[0]].get_query_bed()
    loci = Loci()
    loci.set_use_direction(True)
    loci.set_minimum_distance(args.maximum_intron)
    for i in bestpath:
      r = mpa.entries[i].get_target_bed()
      locus = Locus()
      locus.set_use_direction(True)
      locus.add_member(r)
      loci.add_locus(locus)
    loci.update_loci()
    if len(bestpath) > 1:
      for i in bestpath[1:]:
        query_span = mpa.entries[i].get_query_bed().merge(query_span)
    report = ''
    report += mpa.entries[bestpath[0]].value('qName')+"\t"
    report += str(len(bestpath))+"\t"
    report += str(len(loci.loci))+"\t"
    report += query_span.get_range_string()+"\t"
    report += ','.join([mpa.entries[i].value('strand') for i in bestpath])+"\t"
    report += ','.join([mpa.entries[i].get_query_bed().get_range_string() for i in bestpath])+"\t"
    report += ','.join([mpa.entries[i].get_target_bed().get_range_string() for i in bestpath])+"\t"
    report += ','.join([str(mpa.entries[i].get_quality()) for i in bestpath])+"\t"
    report += ','.join([str(mpa.entries[i].get_coverage()) for i in bestpath])+"\t"
    report += ','.join([str(x) for x in gapsizes])+"\t"
    report += str(besttotalcov)+"\t"
    report += str(bestscore)+"\t"
    report += str(bestsinglescore)+"\t"
    report += str(','.join(query_target_coverages)+"\t")
    #if args.best_report:
    #  best_report_fh.write(report+"\n")
    #for i in bestpath:
    #  args.output.write(mpa.entries[i].get_line()+"\n")
    return [report, [mpa.entries[i].get_line() for i in bestpath]]
    #print '---'
    #print g.get_roots()
if __name__=="__main__":
  main()
