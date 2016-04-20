#!/usr/bin/python
import sys, argparse
import Bio.Format.BGZF

def main():
  parser = argparse.ArgumentParser(description="A gzip compatible bgzf fromat",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="Use - for STDIN")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('-z','--zip',action='store_true',help="compress the file or stream")
  group.add_argument('-x','--unzip',action='store_true',help="uncompress the archive or stream")
  parser.add_argument('-o','--output',help="output file")
  args = parser.parse_args()
  
  of = sys.stdout
  if args.output: of=open(args.output,'wb')
  inf = sys.stdin
  if args.input == '-':
    if args.unzip:
      inf = sys.stdin
      br = Bio.Format.BGZF.reader(inf)
    else: 
      inf = sys.stdin
      bw = Bio.Format.BGZF.writer(of)
  else: 
    if args.unzip:
      inf = open(args.input,'rb')
      br = Bio.Format.BGZF.reader(inf)
    else: 
      inf = open(args.input,'rb')
      bw = Bio.Format.BGZF.writer(of)

  if args.unzip:
    while True:
      v = br.read(1000000)
      if len(v) == 0: break
      of.write(v)
    inf.close()
    of.close()
  else: # we zip it up 
    while True:
      bytes = inf.read(1000000)
      if len(bytes) == 0: break
      if not inf: break
      bw.write(bytes)
    bw.close()
    of.close()

if __name__=="__main__":
  main()
