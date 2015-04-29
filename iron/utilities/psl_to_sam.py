#!/usr/bin/python
import SamBasics, sys, argparse, FileBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('psl',help="FILENAME of psl file (can be gzipped)")
  parser.add_argument('refgenome',help="FASTA of the reference genome")
  parser.add_argument('--fastq_reads',help="FASTQ of the reads")
  parser.add_argument('--fasta_reads',help="FASTA of the reads")
  parser.add_argument('-o',help="FILENAME to save sam output")
  args = parser.parse_args()
  pscf = SamBasics.PSLtoSAMconversionFactory()
  sys.stderr.write("Creating header from reference fasta\n")
  header = SamBasics.construct_header_from_reference_fasta('/Shared/Au/jason/Reference/UCSC/Human/hg19_GRCh37_feb2009/Genome/genome.fa')
  if args.o:
    of = open(args.o,'w')
    of.write(header)
  else:
    sys.stdout.write(header)
  #pscf.construct_header_from_reference_fasta('test.fa')
  sys.stderr.write("Establishing library of reads\n")
  if args.fastq_reads:
    pscf.set_read_fastq(args.fastq_reads)
  elif args.fasta_reads:
    pscf.set_read_fasta(args.fasta_reads)
  sys.stderr.write("Performing conversion\n")
  gfr = FileBasics.GenericFileReader(args.psl)
  while True:
    line = gfr.readline()
    if not line: break
    samline = pscf.convert_line(line.rstrip())
    if args.o:
      of.write(samline+"\n")
    else:
      sys.stdout.write(samline+"\n")
  if args.o:
    of.close()

main()
