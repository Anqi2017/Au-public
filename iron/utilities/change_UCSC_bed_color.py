#!/usr/bin/python
import argparse, sys, re

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="BEDFILE input or use - for STDIN")
  parser.add_argument('color',choices=['black','blue','green','orange','purple','red'])
  args = parser.parse_args()
  color = '0,0,0'
  if args.color == 'black':
    color = '0,0,0'
  if args.color == 'blue':
    color = '67,162,202'
  elif args.color == 'green':
    color = '49,163,84'
  elif args.color == 'orange':
    color = '254,178,76'
  elif args.color == 'purple':
    color = '136,86,167'
  elif args.color == 'red':
    color = '240,59,32'
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  for line in inf:
    if not re.match('\S+\s+\d+\s+\d+',line): #not a bed format print and continue
      print line.rstrip()
      continue
    f = line.rstrip().split("\t")
    if not re.match('\d+,\d+,\d+',f[8]):
      sys.stderr.write("ERROR unexpected format in field 9\n"+line.rstrip()+"\n")
      sys.exit()
    f[8] = color
    print "\t".join(f)

if __name__=="__main__":
  main()
