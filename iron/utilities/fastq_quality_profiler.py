#!/usr/bin/python
import argparse, sys, gzip
import FASTQBasics
from SequenceBasics import FastqHandleReader
import random

def main():
  parser = argparse.ArgumentParser(description="Tool for describing the quality in a fastq")
  parser.add_argument('input',help="FILENAME input fastq or - for stdin (or a profile)")
  parser.add_argument('--quality_type',choices=['S','P','I','J','L'],help="P (Pacbio Phred+33), S (Sanger Phred+33), I (Illumina 1.3+ Phred+64), J (Illumina 1.5+ Phred+64) or autodetect")
  parser.add_argument('--autodetect_depth',type=int,default=100000,help="How many reads to read in on quality detection (stored in memory)")
  parser.add_argument('--training_depth',type=int,help="INT Training read count")
  parser.add_argument('--read_profile',action='store_true')
  args = parser.parse_args()
  inf = sys.stdin
  if args.read_profile: 
    do_reader(args) #jump down to just reading and reporting on the contents of a profile
    return
  if args.input != '-':  
    if args.input[-3:] == '.gz':
      inf = gzip.open(args.input)
    inf = open(args.input)
  fqr = FastqHandleReader(inf)
  buffer = []
  type = None
  if not args.quality_type:
    detector = FASTQBasics.QualityFormatDetector()
    while True:
      entry = fqr.read_entry()
      if not entry: break
      buffer.append(entry)
      detector.record_observation(entry['qual'])
      if len(buffer) >= args.autodetect_depth: break
    type = detector.call_type()
    sys.stderr.write(detector.about+"\n")
  else:
    type = args.quality_type
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
    if z % 100 == 0:
      sys.stderr.write(str(z)+"\r")
  sys.stderr.write("\n")
  print qp.get_serialized()    

def do_reader(args):
  inf = sys.stdin
  if args.input != '-':
    inf = open(args.input)
  line = inf.readline().rstrip()
  inf.close()
  qp = FASTQBasics.QualityProfile()
  qp.read_serialized(line)
  qp.show_stats()
  qconv = FASTQBasics.QualityFormatConverter(qp.quality_type)
  errtot = 0
  tot = 0
  for j in range(2,5):
    rsize = 10**j
    for i in range(0,100):
      qual = qp.emit(rsize)
      for i in range(0,len(qual)):
        p = qconv.call_observed_ascii_probability(qual[i])
        #Cover the special case of 'B' in the illumina 1.5 (J) type
        if qp.quality_type == 'J':
          if qual[i] == 'B':
            nonB = qp.emit_non_B(i,len(qual))
            p = qconv.call_observed_ascii_probability(nonB)
        rnum = random.random()
        if rnum < p: errtot += 1
        tot += 1
    print str(float(errtot)/tot) + "\t" + str(rsize)

if __name__=="__main__":
  main()
