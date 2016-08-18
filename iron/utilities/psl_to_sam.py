#!/usr/bin/python
import SamBasics, sys, argparse, gzip

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('psl',help="FILENAME of psl file (can be gzipped)")
  parser.add_argument('refgenome',help="FASTA of the reference genome")
  parser.add_argument('--min_intron_size',default=68,type=int,help="INT minimum intron size")
  parser.add_argument('--fastq_reads',help="FASTQ of the reads")
  parser.add_argument('--fasta_reads',help="FASTA of the reads")
  parser.add_argument('--skip_directionless_splice',action='store_true',help='only output reads where canonical splice sites indicate direciton if junctions are present')
  parser.add_argument('-o',help="FILENAME to save sam output")
  args = parser.parse_args()
  pscf = SamBasics.PSLtoSAMconversionFactory()
  pscf.set_min_intron_size(args.min_intron_size)
  sys.stderr.write("Creating header from reference fasta\n")
  #header = SamBasics.construct_header_from_reference_fasta('/Shared/Au/jason/Reference/UCSC/Human/hg19_GRCh37_feb2009/Genome/genome.fa')
  if args.skip_directionless_splice:
    pscf.set_skip_directionless_splice()
  header = SamBasics.construct_header_from_reference_fasta(args.refgenome)
  if args.o:
    of = open(args.o,'w')
    of.write(header)
  else:
    sys.stdout.write(header)
  sys.stderr.write("setting reference fasta for conversion\n")
  pscf.set_reference_genome(args.refgenome)
  sys.stderr.write("determining mapping counts from psl\n")
  pscf.set_mapping_counts(args.psl)
  #pscf.construct_header_from_reference_fasta('test.fa')
  sys.stderr.write("Establishing library of reads\n")
  if args.fastq_reads:
    pscf.set_read_fastq(args.fastq_reads)
  elif args.fasta_reads:
    pscf.set_read_fasta(args.fasta_reads)
  sys.stderr.write("Performing conversion\n")
  gfr = None
  if args.psl[-3:]=='.gz': 
    gfr = gzip.open(args.psl)
  else:
    gfr = open(args.psl)
  skipped = 0
  while True:
    line = gfr.readline()
    if not line: break
    samline = pscf.convert_line(line.rstrip())
    if not samline:
      skipped += 1 
      sys.stderr.write("\rskipping directionless splice ("+str(skipped)+")            ")
      continue # happens if we are skipping directionless splice
    if args.o:
      of.write(samline+"\n")
    else:
      sys.stdout.write(samline+"\n")
  if args.o:
    of.close()
  gfr.close()
  sys.stderr.write("\n")
main()
