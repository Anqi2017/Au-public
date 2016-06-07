#!/usr/bin/python
import sys, argparse, gzip, re
from Bio.Format.GPD import GPD
from multiprocessing import cpu_count, Pool, Lock

glock = Lock()
of = sys.stdout
sys.setrecursionlimit(10000)

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  parser.add_argument('--output','-o',help="output file otherwise STDOUT")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Number of threads to convert names")
  parser.add_argument('-r','--reference',required=True,help="reference gpd")
  args = parser.parse_args()
  
  global of
  if args.output:
    if re.search('\.gz$',args.output):
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')

  #read the reference gpd
  rinf = None
  txome = {}
  if re.search('\.gz$',args.reference):
    rinf = gzip.open(args.reference)
  else:
    rinf = open(args.reference)
  sys.stderr.write("Reading in reference\n")
  z = 0
  for line in rinf:
    z += 1
    gpd = GPD(line)
    if z%100 == 0:  sys.stderr.write(str(z)+"          \r")
    if gpd.value('chrom') not in txome: txome[gpd.value('chrom')] = []
    r = gpd.get_range()
    r.set_payload(gpd)
    txome[gpd.value('chrom')].append(r)
  rinf.close()
  sys.stderr.write(str(z)+"          \r")
  sys.stderr.write("\n")
  inf = sys.stdin
  if args.input != '-':
    if re.search('\.gz$',args.input):
      inf = gzip.open(args.input)
    else:
      inf = open(args.input)
  z = 0
  buffer = []
  bsize = 100000
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for line in inf:
    z+=1
    if len(buffer) >= bsize:
      if args.threads > 1:
        p.apply_async(do_buffer,args=(buffer,txome,args),callback=do_out)
      else:
        v = do_buffer(buffer,txome,args)
        do_out(v)
      buffer = []
    buffer.append([line,z])
  if len(buffer) > 0:
      if args.threads > 1:
        p.apply_async(do_buffer,args=(buffer,txome,args),callback=do_out)
      else:
        v = do_buffer(buffer,txome,args)
        do_out(v)
  if args.threads > 1:
    p.close()
    p.join()
  inf.close()
  of.close()

def do_out(results):
  global glock
  global of
  glock.acquire()
  for res in results:
    of.write(res)
  glock.release()

def do_buffer(buffer,txome,args):
  results = []
  for line_z in buffer:
    z = line_z[1]
    line = line_z[0]
    gpd = GPD(line)
    v = annotate_line(gpd,txome,args)
    if not v: continue
    type = 'partial'
    if v[0]: type = 'full'
    exon_count = v[2]    
    most_consecutive_exons = v[3]
    read_exon_count = v[4]
    tx_exon_count = v[5]
    overlap_size = v[6]
    read_length = v[7]
    tx_length = v[8]
    results.append(str(z)+"\t"+gpd.get_gene_name()+"\t"+v[9].get_gene_name()+"\t"+v[9].get_transcript_name()+"\t"+type+"\t"+\
          str(exon_count)+"\t"+str(most_consecutive_exons)+"\t"+str(read_exon_count)+"\t"+str(tx_exon_count)+"\t"+\
          str(overlap_size)+"\t"+str(read_length)+"\t"+str(tx_length)+"\n")
  return results

def annotate_line(gpd,txome,args):
  v = gpd.get_range()
  if v.chr not in txome: return None
  possible = [x.get_payload() for x in txome[v.chr] if x.overlaps(v)]
  candidates = []
  if len(possible) == 0: return None
  for tx in possible:
    eo = None
    full = False
    subset = False
    econsec = 1
    if tx.get_exon_count() == 1 or gpd.get_exon_count() == 1:
      eo = gpd.exon_overlap(tx,single_minover=100,single_frac=0.5)
    else:
      eo = gpd.exon_overlap(tx,multi_minover=10,multi_endfrac=0,multi_midfrac=0.8,multi_consec=False)
      if eo.is_full_overlap():
        full = True
      if eo.is_subset():
        subset = True
      if eo:
        econsec = eo.consecutive_exon_count()
    if not eo: continue
    ecnt = eo.match_exon_count()
    osize = gpd.overlap_size(tx)
    candidates.append([full,subset,ecnt,econsec,gpd.get_exon_count(),tx.get_exon_count(),osize,gpd.get_length(),tx.get_length(),tx])
  if len(candidates)==0: return None
  bests = sorted(candidates,key=lambda x: (-x[0],-x[1],-x[3],-x[2],-min(float(x[6])/float(x[7]),float(x[6])/float(x[8]))))
  return bests[0]

if __name__=="__main__":
  main()
