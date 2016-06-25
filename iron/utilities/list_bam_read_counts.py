#!/usr/bin/python
import sys, argparse
from subprocess import Popen, PIPE
from multiprocessing import cpu_count, Pool

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inputs',nargs='+',help="Use FILE(s)")
  args = parser.parse_args()
  
  p = Pool(processes=cpu_count())
  results = []
  for f in args.inputs:
    v = p.apply_async(check_file,args=(f,))
    results.append(v)
  p.close()
  p.join()
  for res in [x.get() for x in results]:
    print res
def check_file(f):
    cmd1 = 'samtools view '+f
    cmd2 = 'cut -f 1'
    cmd3 = 'sort'
    cmd4 = 'uniq'
    cmd5 = 'wc -l'
    p5 = Popen(cmd5.split(),stdin=PIPE,stdout=PIPE)
    p4 = Popen(cmd4.split(),stdout=p5.stdin,stdin=PIPE)
    p3 = Popen(cmd3.split(),stdout=p4.stdin,stdin=PIPE)
    p2 = Popen(cmd2.split(),stdout=p3.stdin,stdin=PIPE)
    p1 = Popen(cmd1.split(),stdout=p2.stdin)
    p1.communicate()
    p2.communicate()
    p3.communicate()
    p4.communicate()
    return f+"\t"+p5.communicate()[0].rstrip()
    

if __name__=="__main__":
  main()
