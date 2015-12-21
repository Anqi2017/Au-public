#!/usr/bin/python
import argparse, sys
import GTFBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('stringtie_gtf')
  args = parser.parse_args()

  inf = sys.stdin
  if args.stringtie_gtf != '-':
    inf = open(args.stringtie_gtf)
  genes = {}
  for line in inf:
    if line[0] == '#': continue
    gtf = GTFBasics.GTF(line)
    if gtf.entry['gff'][2] != 'transcript': continue
    id = gtf.entry['attributes']['reference_id']
    ref_gene = gtf.entry['attributes']['ref_gene_id']
    gene = gtf.entry['attributes']['gene_id']
    tpm = float(gtf.entry['attributes']['TPM'])
    #print id + "\t" + ref_gene + "\t" + gene + "\t" + str(tpm)
    if ref_gene not in genes:  genes[ref_gene] = {}
    genes[ref_gene][id] = {}
    genes[ref_gene][id]['TPM'] = tpm
    genes[ref_gene][id]['id'] = id
    genes[ref_gene][id]['gene_id'] = gene
    genes[ref_gene][id]['ref_gene_id'] = ref_gene

  for gene in genes:
    tot = 0
    for id in genes[gene]:
      tot += genes[gene][id]['TPM']
    for id in genes[gene]:
      print id + "\t" + gene + "\t" + str(genes[gene][id]['TPM']) + "\t" + str(tot)
    

if __name__=="__main__":
  main()
