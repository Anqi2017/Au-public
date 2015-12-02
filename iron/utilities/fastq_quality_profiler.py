#!/usr/bin/python
import argparse, sys, gzip
import FASTQBasics
from SequenceBasics import FastqHandleReader

def main():
  parser = argparse.ArgumentParser(description="Tool for describing the quality in a fastq")
  parser.add_argument('input',help="FILENAME input fastq or - for stdin")
  parser.add_argument('--quality_type',choices=['S','P','I','J','L'],help="P (Pacbio Phred+33), S (Sanger Phred+33), I (Illumina 1.3+ Phred+64), J (Illumina 1.5+ Phred+64) or autodetect")
  parser.add_argument('--autodetect_depth',type=int,default=100000,help="How many reads to read in on quality detection (stored in memory)")
  parser.add_argument('--training_depth',type=int,help="INT Training read count")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-':  
    if args.input[-3:] == '.gz':
      inf = gzip.open(args.input)
    inf = open(args.input)
  fqr = FastqHandleReader(inf)
  detector = FASTQBasics.QualityFormatDetector()
  buffer = []
  while True:
    entry = fqr.read_entry()
    if not entry: break
    buffer.append(entry)
    detector.record_observation(entry['qual'])
    if len(buffer) >= args.autodetect_depth: break
  type = detector.call_type()
  sys.stderr.write(detector.about+"\n")
  qp = FASTQBasics.QualityProfile(type)
  # Now that we have a type, we can do some profiling
  stats = {}
  z = 0
  while True:
    if args.training_depth:
      if z > args.training_depth: break
    z += 1
    if len(buffer) > 0:
      entry = buffer.pop()
    else:
      entry = fqr.read_entry()
    if not entry: break
    qp.record_observation(entry['qual'])
    if z % 10 == 0:
      sys.stderr.write(str(z)+"\r")
  sys.stderr.write("\n")
  print qp.get_serialized()    

if __name__=="__main__":
  main()
