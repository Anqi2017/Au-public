#!/usr/bin/python 
import argparse, re

def main():
  parser = argparse.ArgumentParser(description="Filter a PSL file by best matches or a name ist",\
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i','--input',required=True,help="FILENAME input")
  parser.add_argument('-o','--output',default='-',help="FILENAME output use '-' for STDOUT")
  parser.add_argument('--names',help="FILENAME a list of names to filter on")
  #parser.add_argument('-v','--invert',help="Reverse the filter to exclude matches")
  parser.add_argument('--best',action='store_true',help="Filter the best entry only")
  args = parser.parse_args()

  nameset = set()
  # read name list if we have one
  if args.names:
    with open(args.names) as inf:
      for name in inf:
        name = name.rstrip()
        nameset.add(name)        
  comments = set()
  keepers = set()
  bestmat = {} # the best matches
  bestz = {} # the line numbers of the best
  z = 0
  with open(args.input) as inf:
    for line in inf:
      line = line.rstrip()
      z += 1
      if re.match('^#',line):  
        comments.add(z)
        continue
      f = line.split('\t')
      matches = int(f[0])
      name = f[9]
      if name not in bestmat:
        bestmat[name] = 0
      if matches > bestmat[name]:
        bestmat[name] = matches
        bestz[name] = z
      if args.names and name in nameset:
        keepers.add(z)
  best_entries = set()
  for name in bestz:
    best_entries.add(bestz[name])
  # now go through a second time
  z = 0
  with open(args.input) as inf:
    for line in inf:
      line = line.rstrip()
      z += 1
      if z in comments:
        print line
        continue
      if args.names and args.best:
        if z in args.best_entries and z in keepers:
          print line
          continue
      if args.names and z in keepers:
        print line
        continue
      if args.best and z in best_entries:
        print line
        continue

if __name__=="__main__":
  main()
