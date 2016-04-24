#!/usr/bin/python
import argparse, sys, re

def main():
  parser = argparse.ArgumentParser(description="output a gene list based on some filtering criteria")
  parser.add_argument('--min_fpkm',type=float,default=0,help="FLOAT minimum FPKM to report")
  parser.add_argument('--length_filter',help="INT,INT i.e. 0,1000 where a sequence of size zero to 1000 will be accetable.")
  parser.add_argument('genes_fpkm_file')
  parser.add_argument('isoforms_fpkm_file')
  args = parser.parse_args()
  genes = read_fpkm(args.genes_fpkm_file)
  isoforms = read_fpkm(args.isoforms_fpkm_file)
  # find any genes with a minimu fpkm meeting our threshold
  expressed_genes = set()
  for tracking_id in genes:
    meets_thresh = 0
    for entry in genes[tracking_id]:
      if entry['fpkm'] >= args.min_fpkm:
        meets_thresh = 1
        break
    if meets_thresh == 0: # this is a bad tracking_id go onto next
      continue
    expressed_genes.add(tracking_id)
  # now go though the isoforms and grab the best length from 
  # each gene and any otheres that meet our thresholds
  best_gene = {}
  for tracking_id in isoforms:
    for entry in isoforms[tracking_id]:
      if entry['gene_name'] not in expressed_genes: continue # don't bother if the gene is not expressed
      if entry['gene_name'] not in best_gene:
        best_gene[entry['gene_name']] = {}
        best_gene[entry['gene_name']]['length'] = 0
        best_gene[entry['gene_name']]['fpkm'] = 0
      if entry['fpkm'] > best_gene[entry['gene_name']]['fpkm']:
        best_gene[entry['gene_name']]['fpkm'] = entry['fpkm']
        best_gene[entry['gene_name']]['length'] = entry['length']
  gene_lengths = {}
  for gene in best_gene:
    if gene not in gene_lengths:
      gene_lengths[gene] = set()
    gene_lengths[gene].add(best_gene[gene]['length'])
  for tracking_id in isoforms:
    for entry in isoforms[tracking_id]:
      if entry['gene_name'] not in expressed_genes: continue # don't bother if the gene is not expressed
      if entry['fpkm'] >= args.min_fpkm: 
        gene_lengths[entry['gene_name']].add(entry['length'])      
  if args.length_filter:
    rvals = [int(x) for x in args.length_filter.split(",")]
    sys.stderr.write("length filter "+str(rvals)+"\n")
    gene_names = gene_lengths.keys()
    for gene in gene_names:
      good_length = 0
      for val in gene_lengths[gene]:
        if val < rvals[0] or val > rvals[1]: continue
        good_length = 1
        break
      if good_length == 0: # we are not in the right length
        del gene_lengths[gene]
  gene_names  = gene_lengths.keys()
  for gene in sorted(gene_names):
    print gene + "\t" + ",".join([str(x) for x in sorted(gene_lengths[gene])])

def read_fpkm(filename):
  entries = {}
  with open(filename) as inf:
    header = inf.readline()
    if not re.match('^tracking_id',header): 
      sys.stderr.write("ERROR expecting header in fpkm file\n")
      sys.exit()
    for line in inf:
      f = line.rstrip().split("\t")
      tracking_id = f[0]
      fpkm = float(f[9])
      gene_name = f[4]
      length = False
      if f[7] != '-':
        length = int(f[7])
      if not tracking_id in entries:
        entries[tracking_id] = []
      entry = {}
      entry['fpkm'] = fpkm
      entry['gene_name'] = gene_name
      entry['length'] = length
      entries[tracking_id].append(entry)
  return entries
main()
