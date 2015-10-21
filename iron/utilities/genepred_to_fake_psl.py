#!/usr/bin/python
import argparse, sys
import GenePredBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input_file',help="use - for STDIN")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input_file != '-':
    inf = open(args.input_file)
  for line in inf:
    e = GenePredBasics.line_to_entry(line.rstrip())
    matches = 0
    qstartslist = []
    for i in range(0,len(e['exonStarts'])):
      mylen = e['exonEnds'][i]-e['exonStarts'][i]
      matches += mylen
      qstartslist.append(matches-mylen)
    qstarts = ','.join([str(x) for x in qstartslist])+','
    oline =  str(matches)+"\t" # 1
    oline += "0\t" # 2
    oline += "0\t" # 3
    oline += "0\t" # 4
    oline += "0\t" # 5
    oline += "0\t" # 6
    oline += "0\t" # 7
    oline += "0\t" # 8
    oline += e['strand']+"\t" # 9
    oline += e['name']+"\t" # 10
    oline += str(matches)+"\t" # 11
    oline += "0\t" # 12
    oline += str(matches)+"\t" # 13
    oline += str(e['chrom'])+"\t" # 14
    oline += str(e['exonEnds'][-1])+"\t" # 15
    oline += str(e['exonStarts'][0])+"\t" # 16
    oline += str(e['exonEnds'][-1])+"\t" # 17
    oline += str(len(e['exonStarts']))+"\t" # 18
    oline += ','.join([str(e['exonEnds'][x]-e['exonStarts'][x]) for x in range(0,len(e['exonStarts']))])+','+"\t" # 19
    oline += qstarts + "\t" # 20
    oline += ','.join([str(x) for x in e['exonStarts']])+',' # 21
    print oline
  inf.close()

if __name__=="__main__":
  main()
