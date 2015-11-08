#!/usr/bin/python
import PSLBasics
import argparse, sys, random, subprocess

def main():
  parser = argparse.ArgumentParser(description="Take a PSL file and output only one alignment per query.  This doesn't do any magic to combine entries, it only seeks the best quality alignment (by most aligned bases).")
  parser.add_argument('input',help="PSLFILE or - for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set.")
  parser.add_argument('--already_ordered',action='store_true',help="The PSL file has already been ordered according to query.  This speeds reading.")
  parser.add_argument('-T','--tempdir',help="DIRECTORY where we can write temporary files.")
  parser.add_argument('-S','--maxtempsize',help="The maximum size for chunks of the temporary files in a sort.  Default units is kb just like linux sort -S... since that is what this feeds the argument to.")
  args = parser.parse_args()
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  # inf2 is our sorted file (by query)
  inf2 = inf
  if not args.already_ordered:
    cmd = "sort -k 10,10"
    if args.tempdir:
      cmd += " -T "+args.tempdir.rstrip('/')
    if args.maxtempsize:
      cmd += " -S "+args.maxtempsize
    p1 = subprocess.Popen(cmd,shell=True,stdin=inf,stdout=subprocess.PIPE)    
    inf2 = p1.stdout
  gr = PSLBasics.GenericOrderedMultipleAlignmentPSLReader()
  gr.set_handle(inf2)
  while True:
    mpa = gr.read_next()
    if not mpa: break
    maxcov = 0
    bestpsl = None
    for psl in mpa.entries:
      cov = psl.get_coverage()
      if cov >= maxcov:
        bestpsl = psl
        maxcov = cov
    of.write(psl.get_line()+"\n")
  of.close()
  if not args.already_ordered:
    p1.communicate() # make sure this is all done
  inf.close()
  
if __name__=="__main__":
  main()
