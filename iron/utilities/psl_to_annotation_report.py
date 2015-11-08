#!/usr/bin/python
import sys, argparse, os, subprocess, re

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
  statof.write(str(query_seqs)+" Aligned reads\n")
  statof.write(str(aligned_query_bases)+" of "+str(query_bases)+" bases in the aligned reads were part of the alignments\n")
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
  statof.write(str(totalreads) + " of "+str(query_seqs)+" reads were annotated\n")
  statof.write(str(multiexon) + " of "+str(totalreads)+" annotated reads were multiexon\n")
  statof.write(str(totalexons)+" of "+str(queryexons)+" exons were annotated from among all aligned exons\n")
  statof.write(str(len(genes))+" genes were identified\n")
  statof.write(str(len(multiexongenes))+" multi-exon genes were identified\n")
  statof.write(str(len(transcripts))+" transcripts were identified\n")
  statof.write(str(len(multiexontranscripts))+" multi-exon transcripts were identified\n")
  statof.close()
  return

if __name__=="__main__":
  main()
