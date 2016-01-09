#!/usr/bin/python 
import argparse, re, sys

def main():
  parser = argparse.ArgumentParser(description="Filter a PSL file by best matches or a name ist",\
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="FILENAME input or - for STDIN")
  parser.add_argument('-o','--output',help="FILENAME output or dont set for STDOUT")
  parser.add_argument('--names',help="FILENAME a list of names to filter on")
  #parser.add_argument('-v','--invert',help="Reverse the filter to exclude matches")
  parser.add_argument('--best',action='store_true',help="Filter the best entry only")
  parser.add_argument('--invert',action='store_true',help="If set, remove these names rather than keep them")
  args = parser.parse_args()

  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')
  nameset = set()
  # read name list if we have one
  if args.names:
    with open(args.names) as inf2:
      for name in inf2:
        name = name.rstrip()
        nameset.add(name)        
  comments = set()
  keepers = set()
  bestmat = {} # the best matches
  bestz = {} # the line numbers of the best
  buffer = []
  z = 0
  for line in inf:
    line = line.rstrip()
    if args.input == '-': buffer.append(line)
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
    if args.names:
      if args.invert:
        if name not in nameset:
          keepers.add(z)
      elif name in nameset:
        keepers.add(z)
  inf.close()
  best_entries = set()
  for name in bestz:
    best_entries.add(bestz[name])
  # now go through a second time
  z = 0
  if args.input == '-': inf = buffer
  else: inf = open(args.input)
  for line in inf:
      line = line.rstrip()
      z += 1
      if z in comments:
        of.write(line+"\n")
        continue
      if args.names and args.best:
        if z in args.best_entries and z in keepers:
          of.write(line+"\n")
          continue
      if args.names and z in keepers:
        of.write(line+"\n")
        continue
      if args.best and z in best_entries:
        of.write(line+"\n")
        continue
  of.close()

if __name__=="__main__":
  main()
