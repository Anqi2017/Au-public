#!/usr/bin/python
import argparse, sys
import GTFBasics

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('stringtie_gtf',help="STRINGTIE GTF FILE or '-' for STDIN")
  parser.add_argument('--by_stringtie_id',action='store_true')
  parser.add_argument('--transcript_coverage',action='store_true')
  args = parser.parse_args()

  inf = sys.stdin
  if args.stringtie_gtf != '-':
    inf = open(args.stringtie_gtf)
  genes = {}
  for line in inf:
    if line[0] == '#': continue
    gtf = GTFBasics.GTF(line)
    if gtf.entry['gff'][2] != 'transcript': continue
    id = None
    if 'reference_id' in gtf.entry['attributes']:
      id = gtf.entry['attributes']['reference_id']
    transcript_id = gtf.entry['attributes']['transcript_id']
    ref_gene = None
    if 'gene_id' in gtf.entry['attributes']:
      ref_gene = gtf.entry['attributes']['gene_id']
    elif 'ref_gene_id' in gtf.entry['attributes']:
      ref_gene = gtf.entry['attributes']['ref_gene_id']
    elif 'ref_gene_name' in gtf.entry['attributes']:
      ref_gene = gtf.entry['attributes']['ref_gene_name']
    gene = gtf.entry['attributes']['gene_id']
    my_gene = ref_gene
    my_id = transcript_id
    if args.by_stringtie_id: 
      my_gene = gene
      my_id = transcript_id
    if not args.by_stringtie_id and not my_gene:
      sys.stderr.write("ERROR: this is a predicted set, you cannot use only a reference annotation for it.  try --by_stringtie_id\n"+line+"\n")
      sys.exit()
    tpm = float(gtf.entry['attributes']['TPM'])
    cov = float(gtf.entry['attributes']['cov'])
    #print id + "\t" + ref_gene + "\t" + gene + "\t" + str(tpm)
    if my_gene not in genes:  genes[my_gene] = {}
    genes[my_gene][my_id] = {}
    genes[my_gene][my_id]['TPM'] = tpm
    genes[my_gene][my_id]['cov'] = cov
    genes[my_gene][my_id]['ref_id'] = id
    genes[my_gene][my_id]['transcript_id'] = transcript_id
    genes[my_gene][my_id]['gene_id'] = gene
    genes[my_gene][my_id]['ref_gene_id'] = my_gene
  if args.transcript_coverage:
    for gene in genes:
      for transcript in genes[gene]:
        print transcript+"\t"+str(genes[gene][transcript]['cov'])
    return

  for gene in genes:
    tot = 0
    for id in genes[gene]:
      tot += genes[gene][id]['TPM']
    for id in genes[gene]:
      if args.by_stringtie_id:
        ref_tx = ''
        ref_gene = ''
        if genes[gene][id]['ref_id']: ref_tx =  genes[gene][id]['ref_id']
        if genes[gene][id]['ref_gene_id']: ref_gene =  genes[gene][id]['ref_gene_id']
        print genes[gene][id]['gene_id'] + "\t" + genes[gene][id]['transcript_id'] + "\t" + str(genes[gene][id]['TPM']) + "\t" + str(tot) + "\t" + ref_tx + "\t" + ref_gene
      else:  
        print id + "\t" + gene + "\t" + str(genes[gene][id]['TPM']) + "\t" + str(tot)

if __name__=="__main__":
  main()
