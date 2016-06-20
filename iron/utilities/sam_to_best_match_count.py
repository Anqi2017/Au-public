#!/usr/bin/python
import sys, argparse
from SamBasics import SamStream, SAM
from subprocess import Popen, PIPE
from multiprocessing import Pool, Lock

glock = Lock()
best_matches = {}

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="BAM Use - for sam STDIN")
  parser.add_argument('--threads',type=int,default=1,help="Number of threads, use 1 if input is stdin")
  args = parser.parse_args()
  multi = None
  if args.input == '-':
    args.input = [SamStream(sys.stdin)]
    args.threads = 1
  else:
    inputs = []
    multi = []
    for i in range(0,args.threads):
      cmd = 'samtools view -h '+args.input
      p = Popen(cmd.split(),stdout=PIPE)
      inputs.append(SamStream(p.stdout))
      multi.append(p)
    args.input = inputs

  if args.threads > 1:
    q = Pool(processes=args.threads)
  for i in range(0,args.threads):
    if args.threads > 1:
      q.apply_async(get_best,args=(args.input[i],i,args.threads),callback=mergebest)
    else:
      v = get_best(args.input[i],i,args.threads)
      mergebest(v)
  if args.threads > 1:
    q.close()
    q.join()

  #now can output results
  global best_matches
  for name in best_matches:
    print name +"\t"+str(best_matches[name])

  if multi:
    for p in multi:
      p.communicate()

def mergebest(v):
  global best_matches
  global glock
  glock.acquire()
  for name in v:
    if name not in best_matches:
      best_matches[name] = 0
    if v[name] > best_matches[name]:
      best_matches[name] = v[name]
  glock.release()
  return

def get_best(input,i,total):
  best = {}
  cnt = 0
  while True:
    if cnt%total != i: continue
    s = input.read_entry()
    if not s: break
    s = SAM(s)
    if s.entry['qname'] not in best:
      best[s.entry['qname']] = s.get_coverage()
    elif s.get_coverage() > best[s.entry['qname']]:
      best[s.entry['qname']] = s.get_coverage()
  return best

if __name__=="__main__":
  main()
