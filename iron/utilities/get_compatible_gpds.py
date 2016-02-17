#!/usr/bin/python
import sys, argparse
from GenePredBasics import GenePredEntry as GPD
from GenePredFuzzyBasics import FuzzyGenePred
from multiprocessing import Pool, cpu_count, Lock

glock = Lock()
of = None

def main():
  parser = argparse.ArgumentParser(description="do join of gpd overlaps")
  parser.add_argument('gpd_left')
  parser.add_argument('gpd_right')
  parser.add_argument('-j','--junction_tolerance',default=5,type=int)
  parser.add_argument('--threads',type=int,default=1)
  parser.add_argument('--left_outer_join',action='store_true')
  parser.add_argument('-o','--output')
  args = parser.parse_args()

  global of
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  
  #read the left into memory
  rightlines = []
  z = 0
  with open(args.gpd_right) as inf:
    for line in inf:
      z += 1
      if z%100==0: sys.stderr.write("reading in "+str(z)+"    \r")
      fgpd = FuzzyGenePred(GPD(line.rstrip()),juntol=args.junction_tolerance)
      rightlines.append(fgpd)
  sys.stderr.write("\n")
  sys.stderr.write("finished reading in "+str(len(rightlines))+" gpd entries\n")
  i = 0
  if args.threads > 1:
    p = Pool(processes=args.threads)
  with open(args.gpd_left) as inf:
    for line in inf:
      i += 1
      fgpd1 = FuzzyGenePred(GPD(line.rstrip()),juntol=args.junction_tolerance)
      if args.threads > 1:
        p.apply_async(get_compatible,args=(rightlines,fgpd1,i,args),callback=do_output)
      else:
        v = get_compatible(rightlines,fgpd1,i,args)
        do_output(v)
  if args.threads > 1:
    p.close()
    p.join()
  of.close()
def do_output(v):
  if not v: return
  global glock
  global of
  glock.acquire()
  sys.stderr.write("processing "+v.rstrip().split("\t")[0]+"   \r")
  of.write(v.rstrip()+"\n")
  glock.release()

def get_compatible(rightlines,fgpd1,i,args):
  output =  ''
  j = 0
  cnt = 0
  for fgpd2 in rightlines:
    j += 1
    if len(fgpd1.fuzzy_junctions) != len(fgpd2.fuzzy_junctions): continue
    if fgpd2.compatible_overlap(fgpd1) and len(fgpd1.fuzzy_junctions) == len(fgpd2.fuzzy_junctions):
      # should be exact overlaps
      output += str(i) + "\t" + fgpd1.gpds[0].value('gene_name') + "\t" + fgpd1.gpds[0].value('name') + "\t" + str(j)+"\t" + fgpd2.gpds[0].value('gene_name') + "\t" + fgpd2.gpds[0].value('name')+"\n"
      cnt += 1
  if cnt == 0 and args.left_outer_join: # its an outer join
    output += str(i)+"\t"+fgpd1.gpds[0].value('gene_name')+"\t"+fgpd1.gpds[0].value('name')+"\t"+"0"+"\t"+"\t"
  if output == '': return None
  return output

if __name__=="__main__":
  main()
