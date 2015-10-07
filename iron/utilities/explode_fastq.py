#!/usr/bin/python
import os, sys, argparse
from FileBasics import GenericFileReader

def main():
  parser = argparse.ArgumentParser(description='Split FASTQ file(s) into smaller ones with as many entries as you specify')
  parser.add_argument('size',type=int,help='Number of sequences to put into each file')
  parser.add_argument('output_directory',help='Name of the directory to put sequences')
  parser.add_argument('fastq_files',nargs='+',help='FILENAME(S) for fastq files')
  args = parser.parse_args()
  if len(args.fastq_files) > 2:
    sys.stderr.write("ERROR only two fastq files at most are supported\n")
    return
  if os.path.exists(args.output_directory):
    sys.stderr.write("ERROR output directory exists already\n")
    return
  os.makedirs(args.output_directory)
  if len(args.fastq_files) == 1:
    out_iter = 1
    fcount = 0
    of = open(args.output_directory.rstrip('/')+'/'+str(out_iter)+'.fq','w')
    gfr = sys.stdin
    if(args.fastq_files[0] != '-'):
      gfr = GenericFileReader(args.fastq_files[0])
    while True:
      lineA = gfr.readline()
      if not lineA: break
      lineB = gfr.readline()
      lineC = gfr.readline()
      lineD = gfr.readline()
      of.write(lineA)
      of.write(lineB)
      of.write(lineC)
      of.write(lineD)
      fcount += 1
      if args.size <= fcount:
        fcount = 0
        out_iter += 1
        of.close()
        of = open(args.output_directory.rstrip('/')+'/'+str(out_iter)+'.fq','w')   
    gfr.close()
  else: # we have two fastq files
    out_iter = 1
    fcount = 0
    of1 = open(args.output_directory.rstrip('/')+'/'+str(out_iter)+'_1.fq','w')
    gfr1 = GenericFileReader(args.fastq_files[0])
    of2 = open(args.output_directory.rstrip('/')+'/'+str(out_iter)+'_2.fq','w')
    gfr2 = GenericFileReader(args.fastq_files[1])
    while True:
      line1a = gfr1.readline()
      line2a = gfr2.readline()
      if not line1a or not line2a: 
        if line1a or line2a:
          sys.stderr.write("WARNING paired file lengths appear different\n")
        break
      line1b = gfr1.readline()
      line2b = gfr2.readline()
      line1c = gfr1.readline()
      line2c = gfr2.readline()
      line1d = gfr1.readline()
      line2d = gfr2.readline()
      of1.write(line1a)
      of2.write(line2a)
      of1.write(line1b)
      of2.write(line2b)
      of1.write(line1c)
      of2.write(line2c)
      of1.write(line1d)
      of2.write(line2d)
      fcount += 1
      if args.size <= fcount:
        fcount = 0
        out_iter += 1
        of1.close()
        of2.close()
        of1 = open(args.output_directory.rstrip('/')+'/'+str(out_iter)+'_1.fq','w')   
        of2 = open(args.output_directory.rstrip('/')+'/'+str(out_iter)+'_2.fq','w')   
    gfr1.close()
    gfr2.close()

main()
