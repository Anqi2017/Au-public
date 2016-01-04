#!/usr/bin/python
import argparse, sys, random, os, subprocess
import SequenceBasics
from shutil import rmtree

# Wrapper for a pipeline to convert BWA-mem results into a psl file


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("bwa_bam",help="BAMFILE or - for sam streamed to stdin")
  parser.add_argument('-o','--output',help="OUTFILE or if not set use STDOUT")
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument("--query_fasta",help="FASTA for query sequences")
  group1.add_argument("--query_fastq",help="FASTQ for query sequences")
  parser.add_argument("--ref",required=True,help="FASTA for reference genome")
  parser.add_argument("-S","--size",help="linux sort option S, in kb if units not specified")
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument("--tempdir",default='/tmp',help="DIR to store a temp directory")
  group2.add_argument("--specific_tempdir",help="Exact DIR to work in, is not deleted")
  parser.add_argument("--get_all_alignments",action='store_true',help="Be exhaustive in retrieving secondary and chimeric alignments")
  args = parser.parse_args()

  # Manage your tempdir
  # put it into args.tempdir
  if not args.specific_tempdir:
    rnum = random.randint(1,100000000)
    args.tempdir = args.tempdir.rstrip('/')+'/weirathe.'+str(rnum) 
  else:
    args.tempdir = args.specific_tempdir.rstrip('/')
  if not os.path.exists(args.tempdir):
    os.makedirs(args.tempdir)

  inf = sys.stdin
  if args.bwa_bam != '-':
    p1 = subprocess.Popen(("samtools view "+args.bwa_bam).split(),stdout=subprocess.PIPE)
    inf = p1.stdout
  # 1. Now we can start the process.  First convert the sam to a psl
  cmd = "sam_to_psl.py -"
  if args.get_all_alignments:
    cmd += " --get_all_alignments"

  #No matter what we don't need to see the exact same alignment more than once
  cmd += " | sort -T "+args.tempdir 
  if args.size:
    cmd += " -S "+args.size
  cmd += " | uniq"

  #Sort the alignment based on query name
  cmd += " | sort_psl.py - --tempdir "+args.tempdir
  if args.size:
    cmd += " -S "+args.size

  #Put the fragmented alignments together
  cmd += " | defragment_PSL_alignments.py - | sort_psl.py - --tempdir "+args.tempdir+" -o "+args.tempdir+"/1.psl"
  if args.size:
    cmd += " -S "+args.size
  p2 = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE)
  for line in inf:
    p2.stdin.write(line)
  p2.communicate()
  if args.bwa_bam != '-':
    p1.communicate()

  # 2. Get a sorted fasta file
  cmd = "sort_fasta.py - --tempdir "+args.tempdir+" -o "+args.tempdir+'/2.query.fa'
  if args.size:
    cmd += " -S "+args.size
  p3 = subprocess.Popen(cmd.split(),stdin=subprocess.PIPE)
  if args.query_fasta:
    gfr = SequenceBasics.GenericFastaFileReader(args.query_fasta)
  else:
    gfr = SequenceBasics.GenericFastqFileReader(args.query_fastq)
  while True:
    e = gfr.read_entry()
    if not e: break
    p3.stdin.write('>'+e['name']+"\n"+e['seq']+"\n")
  p3.communicate()
  gfr.close()

  # 3. Fix the psl and output the results
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')
  cmd = "fix_psl_stats.py "+args.tempdir+'/1.psl '+args.ref+' '+args.tempdir+'/2.query.fa | sort_psl.py - --tempdir '+args.tempdir
  if args.size:
    cmd += " -S "+args.size
  p4 = subprocess.Popen(cmd,shell=True,stdout=of)
  p4.communicate()

  # Clean up your temporary directory if you aren't in a specific one.
  if not args.specific_tempdir:
    rmtree(args.tempdir)
  return

if __name__=="__main__":
  main()

