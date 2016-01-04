#!/usr/bin/python
import argparse, sys, re
import SamBasics
from SequenceBasics import rc
from PSLBasics import PSL

def main():
  parser = argparse.ArgumentParser(description="Convert a sam file into a psl file")
  parser.add_argument('--genome',help="FASTA input file of reference genome")
  parser.add_argument('--get_secondary_alignments',action='store_true',help="Report SA:Z secondary alignments as well")
  parser.add_argument('--get_alternative_alignments',action='store_true',help="Report XA:Z alternative alignments as well")
  parser.add_argument('--get_all_alignments',action='store_true',help="Report SA:Z and XA:Z alternative alignments as well")
  parser.add_argument('--give_unique_names',action='store_true',help="Output query names will be unique.")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--output_fasta',help="FILENAME to save an outgoing fasta.  Only works for primary alignments.")
  group.add_argument('--output_fastq',help="FILENAME to save an outgoing fastq.  Only works for primary alignments.")
  parser.add_argument('infile',help="FILENAME input file or '-' for STDIN")
  parser.add_argument('-o','--output',help="FILENAME for the output, STDOUT if not set.")
  args = parser.parse_args()
  if (args.output_fasta or args.output_fastq) and (args.get_secondary_alignments or args.get_alternative_alignments or args.get_all_alignments):
    sys.stderr.write("ERROR, can only output the fastq/fasta if we are doing primary alignments only.\n")
    sys.exit()
  inf = sys.stdin
  if args.infile != '-': 
    inf = open(args.infile)
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')
  spcf = SamBasics.SAMtoPSLconversionFactory()
  if args.genome: spcf.set_genome(args.genome)
  off = None
  if args.output_fasta:
    off = open(args.output_fasta,'w')
  if args.output_fastq:
    off = open(args.output_fastq,'w')
  z = 0
  for line in inf:
    line = line.rstrip()
    if SamBasics.is_header(line): 
      spcf.read_header_line(line)
      continue
    # We have a line to convert
    psl = spcf.convert_line(line)
    if psl:
      pobj = PSL(psl)
      z += 1
      if args.give_unique_names:
        pobj.entry['qName'] = 'Q'+str(z)
      of.write(pobj.get_line()+"\n")
      if args.output_fastq or args.output_fasta:
        sam = SamBasics.SAM(line)
        sequence = sam.value('seq').upper()
        quality = sam.value('qual')
        if sam.check_flag(16):
          sequence = rc(sam.value('seq').upper())
          quality = sam.value('qual')[::-1]
        if args.output_fasta:
          off.write(">"+pobj.value('qName')+"\n"+sequence+"\n")
        elif args.output_fastq:
          if len(sequence) == len(quality):
            off.write("@"+pobj.value('qName')+"\n"+sequence+"\n"+"+\n"+quality+"\n")
          else:
            sys.stderr.write("ERROR: sequence "+sequence+" length ("+str(len(sequence))+") doesnt match quality "+quality+" length ("+str(len(quality))+")\n")
            sys.exit()
    # Lets look for secondary alignments to convert
    if args.get_secondary_alignments or args.get_all_alignments:
      secondary_alignments = SamBasics.get_secondary_alignments(line.rstrip())
      for samline in secondary_alignments:
        psl = spcf.convert_line(samline)
        if psl:
          #print "\nsecondary"
          #print samline
          z += 1
          pobj = PSL(psl)
          if args.give_unique_names:
            pobj.entry['qName'] = 'Q'+str(z)
          of.write(pobj.get_line()+"\n")
    if args.get_alternative_alignments or args.get_all_alignments:
      alternative_alignments = SamBasics.get_alternative_alignments(line.rstrip())
      for samline in alternative_alignments:
        psl = spcf.convert_line(samline)
        if psl:
          #print "\nsecondary"
          #print samline
          z += 1
          pobj = PSL(psl)
          if args.give_unique_names:
            pobj.entry['qName'] = 'Q'+str(z)
          of.write(pobj.get_line()+"\n")
  inf.close()
  of.close()


if __name__=="__main__":
  main()
