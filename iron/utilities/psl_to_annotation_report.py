#!/usr/bin/python
import sys, argparse, os, subprocess, re
import SequenceBasics

# Pre: Take a psl file and genepred file 
# Post: Create a folder containing the best psl alignments
#       and the genepred equivlents, and the gene annotations 
#       associated with those alignments in best and raw form

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input_psl')
  parser.add_argument('input_gpd')
  parser.add_argument('output_dir')
  parser.add_argument('--run_name',help="Name of the run or will use the PSL file name")
  parser.add_argument('--minoverlap',help="minoverlap parameter of annotated_psl_with_gpd.py")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--query_fasta')
  group.add_argument('--query_fastq')
  args = parser.parse_args()
  gpd_basename = args.input_gpd
  m = re.search('([^\/]+)$',args.input_gpd)
  if m:  gpd_basename = m.group(1)
  args.output_dir = args.output_dir.rstrip('/')
  if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)
  logof = open(args.output_dir+'/log.txt','w')
  cmd = "psl_to_best_psl_entries.py "+args.input_psl+" -o "+args.output_dir+"/best.psl"
  logof.write(cmd+"\n")
  subprocess.call(cmd.rstrip().split())
  cmd =  "annotate_psl_with_gpd.py "+args.output_dir+'/best.psl '+args.input_gpd
  cmd += ' --bestoutput '+args.output_dir+'/best-annotations.txt'
  cmd += ' --rawoutput '+args.output_dir+'/raw-annotations.txt'
  cmd += ' --output '+args.output_dir+'/report-annotations.txt'
  cmd += ' --gpdout '+args.output_dir+'/best.gpd'
  if args.minoverlap:
    cmd += ' --minoverlap '+args.minoverlap
  logof.write(cmd+"\n")
  subprocess.call(cmd.rstrip().split())
  # Make a file condusive to sql import
  bestof = open(args.output_dir+'/best-annotations.forSql','w')
  run_name = args.input_psl
  m = re.search('([^\/]+)$',args.input_psl)
  if m: run_name = m.group(1)
  if args.run_name:  run_name = args.run_name
  with open(args.output_dir+'/best-annotations.txt') as inf:
    for line in inf:
      if re.match('^psl_entry_id',line): continue # skip the header
      bestof.write(run_name+"\t"+line.rstrip()+"\n")
  bestof.close()
  cmd = "genepred_to_bed.py "+args.output_dir+"/best.gpd --headername "+run_name+" > "+args.output_dir+"/for-UCSC.bed"
  logof.write(cmd+"\n")
  subprocess.call(cmd,shell=True)
  if not os.path.exists(args.output_dir+'/gene_rarefraction'):
    os.makedirs(args.output_dir+'/gene_rarefraction')
  cmd = "bestannotation_into_rarefraction.py "+args.output_dir+"/best-annotations.txt --annotation_type gene --gpd_name "+gpd_basename+" > "+args.output_dir+"/gene_rarefraction/rarefraction-genes-all.txt"
  logof.write(cmd+"\n")
  subprocess.call(cmd,shell=True)
  if not os.path.exists(args.output_dir+'/transcript_rarefraction'):
    os.makedirs(args.output_dir+'/transcript_rarefraction')
  cmd = "bestannotation_into_rarefraction.py "+args.output_dir+"/best-annotations.txt --annotation_type gene --gpd_name "+gpd_basename+" > "+args.output_dir+"/transcript_rarefraction/rarefraction-transcripts-all.txt"
  logof.write(cmd+"\n")
  subprocess.call(cmd,shell=True)
  logof.close()
  # Read in the query's if they were specified
  original_queries = None
  original_bp = None
  if args.query_fastq or args.query_fasta:
    original_queries = 0
    original_bp = 0
    gfr = None
    if args.query_fasta:
      gfr = SequenceBasics.GenericFastaFileReader(args.query_fasta)
    elif args.query_fastq:
      gfr = SequenceBasics.GenericFastqFileReader(args.query_fastq)
    while True:
      e = gfr.read_entry()
      if not e: break
      original_queries += 1
      original_bp += len(e['seq'])
    gfr.close()
  statof = open(args.output_dir+'/stats.txt','w')
  # now collect information from these files to make the stats report.
  #1. Number of query sequences aligned
  #2. Number of aligned query bases \t Number of query bases aligned
  #3. Number of multi-exon aligned sequences
  #4. Total number of exons
  #5. Number of genes identified \t Number of genes identifed by multi-exon alignments
  #6. Number of transcripts identified \t Number of transcripts identifed by multi-exon alignments
  query_seqs = 0
  query_bases = 0
  aligned_query_bases = 0
  with open(args.output_dir+'/best.psl') as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      query_seqs += 1
      query_bases += int(f[10])
      aligned_query_bases += int(f[0])+int(f[1]) # matches and mismatches
  statof.write(percent(query_seqs,original_queries)+" " +str(query_seqs)+" Aligned reads")
  if original_queries:
    statof.write(" of "+str(original_queries)+" Total reads")
  statof.write("\n")
  if original_bp:
    statof.write(percent(aligned_query_bases,original_bp)+" "+str(aligned_query_bases)+" of "+str(original_bp)+" aligned bases out of the total bases\n")
  statof.write(percent(aligned_query_bases,query_bases)+" "+str(aligned_query_bases)+" of "+str(query_bases)+" aligned bases in the aligned reads\n")
  queryexons = 0
  with open(args.output_dir+'/best.gpd') as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      queryexons += int(f[8])
  totalexons = 0
  full = 0
  totalreads = 0
  multiexon = 0
  genes = set()
  transcripts = set()
  multiexongenes = set()
  multiexontranscripts = set()
  with open(args.output_dir+'/best-annotations.txt') as inf:
    header = inf.readline()
    for line in inf:
      f = line.rstrip().split("\t")
      totalreads += 1
      genes.add(f[6])
      transcripts.add(f[7])
      if f[9] == 'Full':
        full+=1
      totalexons += int(f[12])
      if int(f[12]) > 1: 
        multiexon += 1
        multiexongenes.add(f[6])
        multiexontranscripts.add(f[7])
  if original_queries:
    statof.write(percent(totalreads,original_queries)+" "+str(totalreads) + " of "+str(original_queries)+" of the total reads were annotated\n")
  statof.write(percent(totalreads,query_seqs)+" "+str(totalreads) + " of "+str(query_seqs)+" of aligned reads were annotated\n")
  statof.write(percent(multiexon,totalreads)+" "+str(multiexon) + " of "+str(totalreads)+" annotated reads were multiexon\n")
  statof.write(percent(totalexons,queryexons)+ " " +str(totalexons)+" of "+str(queryexons)+" exons were annotated from among all aligned exons\n")
  statof.write(str(len(genes))+" genes were identified\n")
  statof.write(str(len(multiexongenes))+" multi-exon genes were identified\n")
  statof.write(str(len(transcripts))+" transcripts were identified\n")
  statof.write(str(len(multiexontranscripts))+" multi-exon transcripts were identified\n")
  statof.close()
  return

def percent(num1,num2):
  num1 = float(num1)
  num2 = float(num2)
  if num2 == 0:
    return "None"
  return str("{0:.3f}%".format(100*num1/num2))

if __name__=="__main__":
  main()
