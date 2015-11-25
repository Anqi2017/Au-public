#!/usr/bin/python
import argparse, sys, re, random, os, gzip
from subprocess import Popen, PIPE
import SimulationBasics
from VCFBasics import VCF
from SequenceBasics import read_fasta_into_hash
from TranscriptomeBasics import Transcriptome

def main():
  parser = argparse.ArgumentParser(description="Create a simulated RNA-seq dataset")
  parser.add_argument('reference_genome',help="The reference genome.")
  parser.add_argument('transcripts_genepred',help="A genepred file describing the transcripts.  Each transcript name must be unique.")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--uniform_expression',action='store_true',help="Uniform distribution of transcript expression")
  group.add_argument('--isoform_expression',help="The transcript expression in TSV format <Transcript name> tab <Expression>")
  group.add_argument('--cufflinks_isoform_expression',help="The expression of the isoforms or - for a uniform distribution of transcript expression")
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--long_reads_only',action='store_true')
  group2.add_argument('--short_reads_only',action='store_true')
  group2.add_argument('--output',help="Directory name for output")
  parser.add_argument('--short_read_count',type=int,default=10000,help="INT number of short reads")
  parser.add_argument('--short_read_length',type=int,default=101,help="INT length of the short reads")
  parser.add_argument('--long_read_count',type=int,default=4000,help="INT default number of long reads")
  args = parser.parse_args()
  if args.output:
    args.output = args.output.rstrip('/')
  ref = read_fasta_into_hash(args.reference_genome)
  txn = Transcriptome()
  txn.set_reference_genome_dictionary(ref)
  with open(args.transcripts_genepred) as inf:
    for line in inf:
      if line[0]=='#': continue
      txn.add_genepred_line(line.rstrip())
  if args.isoform_expression:
    sys.stderr.write("Reading expression from a TSV\n")
    with open(args.isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        txn.add_expression(f[0],float(f[1]))
  elif args.uniform_expression:
    sys.stderr.write("Using uniform expression model\n")
  elif args.cufflinks_isoform_expression:
    sys.stderr.write("Using cufflinks expression\n")
    with open(args.cufflinks_isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        txn.add_expression(f[0],float(f[9]))
  # Now we have the transcriptomes set
  #Now our dataset is set up
  rbe = SimulationBasics.RandomTranscriptomeEmitter(txn)
  if args.short_reads_only:
    rbe.set_gaussian_fragmentation_default_hiseq()
    for zi in range(0,args.short_read_count):
      [name,seq] = rbe.emit_short_read(args.short_read_length)
      print "@SRSIM"+str(zi+1)
      print seq
      print "+"
      print 'I'*len(seq)
    return
  if args.long_reads_only:
    rbe.set_gaussian_fragmentation_default_pacbio()
    for zi in range(0,args.long_read_count):
      [name,seq] = rbe.emit_long_read()
      g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(zi+1)+'/ccs'
      print "@"+g
      print seq
      print "+"
      print 'I'*len(seq)    
    return
  rbe.set_gaussian_fragmentation_default_hiseq()
  # Lets prepare to output now
  if not os.path.exists(args.output):
    os.makedirs(args.output)
  sys.stderr.write("Sequencing short reads\n")
  of1 = gzip.open(args.output+"/SR_1.fq.gz",'wb')
  of2 = gzip.open(args.output+"/SR_2.fq.gz",'wb')
  for i in range(0,args.short_read_count):
    z = i+1
    if z %100==0: sys.stderr.write(str(z)+"\r")
    [name,l1,r1] = rbe.emit_paired_short_read(args.short_read_length)
    of1.write("@SRSIM"+str(z)+"\n")
    of1.write(l1+"\n")
    of1.write("+\n")
    of1.write(len(l1)*'I'+"\n")
    of2.write("@SRSIM"+str(z)+"\n")
    of2.write(r1+"\n")
    of2.write("+\n")
    of2.write(len(r1)*'I'+"\n")
  sys.stderr.write("\nFinished sequencing short reads\n")
  of1.close()
  of2.close()

  # Now lets print out some of the emission details
  of = open(args.output+"/SR_report.txt",'w')
  for name in sorted(rbe.emissions_report.keys()):
    express = 1
    if rbe.transcriptome.expression:
      express = rbe.transcriptome.expression.get_expression(name)
    of.write(name +"\t"+str(express)+"\t"+str(rbe.emissions_report[name])+"\n")
  of.close()
  rbe.emissions_report = {}

  # Now lets create the long read set
  rbe.set_gaussian_fragmentation_default_pacbio()
  sys.stderr.write("Sequencing long reads\n")
  of = gzip.open(args.output+"/LR.fq.gz",'wb')
  for i in range(0,args.long_read_count):
    z = i+1
    if z %100==0: sys.stderr.write(str(z)+"\r")
    [name,seq] = rbe.emit_long_read()
    g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
    of.write("@"+g+"\n")
    of.write(seq+"\n")
    of.write("+\n")
    of.write(len(seq)*'I'+"\n")
  of.close()
  sys.stderr.write("\nFinished sequencing long reads\n")
  # Now lets print out some of the emission details
  of = open(args.output+"/LR_report.txt",'w')
  for name in sorted(rbe.emissions_report.keys()):
    express = 1
    if rbe.transcriptome.expression:
      express = rbe.transcriptome.expression.get_expression(name)
    of.write(name +"\t"+str(express)+"\t"+str(rbe.emissions_report[name])+"\n")
  of.close()
  combo = {}
  with open(args.output+"/SR_report.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      [name,express,left] = f
      if name not in combo:
        combo[name] = {}
        combo[name]['express'] = express
        combo[name]['left'] = 0
      combo[name]['left'] += int(left)
  with open(args.output+"/LR_report.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      [name,express,left] = f
      if name not in combo:
        combo[name] = {}
        combo[name]['express'] = express
        combo[name]['left'] = 0
      combo[name]['left'] += int(left)
  of = open(args.output+"/LR_SR_combo_report.txt",'w')
  for name in sorted(combo):
    of.write(name+"\t"+combo[name]['express']+"\t"+str(combo[name]['left'])+"\n")
  of.close()
if __name__=="__main__":
  main()
