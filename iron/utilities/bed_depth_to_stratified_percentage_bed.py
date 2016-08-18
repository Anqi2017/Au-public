#!/usr/bin/python
import sys, argparse, gzip

g_version = None

def main(args):

  
  inf = sys.stdin
  if args.input != '-': 
    if args.input[-3:]=='.gz': inf = gzip.open(args.input)
    else: inf = open(args.input)

  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz': of = gzip.open(args.output,'w')
    else: of = open(args.output,'w')
  depths = {}
  vals = []
  z = 0
  for line in inf:
    z += 1
    if z % 100000 == 0: sys.stderr.write(str(z)+"    bed entries read   \r")
    f = line.rstrip().split("\t")
    vals.append([f[0],int(f[1]),int(f[2]),int(f[3])])
  z = 0
  sys.stderr.write("\n")
  for f in vals:
    z += 1
    if z % 100000 == 0: sys.stderr.write(str(z)+"    bed entries read   \r")
    #keep track of the number of bases at each depth
    depth = int(f[3])
    cov = int(f[2])-int(f[1])
    if depth not in depths:  depths[depth] = 0
    depths[depth] += cov
    #vals.append([f[0],int(f[1]),int(f[2]),depth])
  sys.stderr.write("\n")
  total_bases = sum(depths.values())
  thresh = {}
  for i in range(1,args.strata+1):
    pos = 0
    cur = float(i)*float(total_bases)/float(args.strata)
    for d in sorted(depths.keys()):
      pos += depths[d]
      if float(pos) > cur:
        thresh[d] = [pos,i]
        break
  cutoffs = sorted(thresh.keys())
  z = 0
  for val in vals:
    z += 1
    if z % 100000 == 0: sys.stderr.write(str(z)+"    bed entries read   \r")
    #print val
    ostrat = args.strata
    for cutoff in cutoffs:
      if val[3] <= cutoff:
        ostrat = thresh[cutoff][1]
        break
    val[3] = ostrat
  sys.stderr.write("\n")
  buffer = vals[0]
  for val in vals[1:]:
    if val[1]==buffer[2] and val[3]==buffer[3] and val[0]==buffer[0]:
      #print 'hello'
      buffer[2] = val[2]
      continue
    else:
      of.write(buffer[0]+"\t"+str(buffer[1])+"\t"+str(buffer[2])+"\t"+str(buffer[3])+"\n")
      buffer = val
  of.write(buffer[0]+"\t"+str(buffer[1])+"\t"+str(buffer[2])+"\t"+str(buffer[3])+"\n")
  of.close()


def do_inputs():
  parser = argparse.ArgumentParser(description="convert bed depth from depth to strata",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="bed or - for STDIN")
  parser.add_argument('strata',type=int,help="number of strata to group reads into")
  parser.add_argument('-o','--output',help="output file")
  args = parser.parse_args()
  return args

def external_cmd(cmd,version=None):
  #set version by input
  global g_version
  g_version = version

  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv
  
if __name__=="__main__":
  args = do_inputs()
  main(args)
