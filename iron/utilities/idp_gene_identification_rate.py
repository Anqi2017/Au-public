#!/usr/bin/python
import argparse, subprocess, os, multiprocessing, re, sys
from random import randint
from GenePredBasics import line_to_entry as genepred_line_to_entry
from shutil import copytree, rmtree

def main():
  parser = argparse.ArgumentParser(description="Generate figure for gene identification rates in IDP.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('short_read_alignment',help="FILENAME short read alignment.")
  parser.add_argument('transcriptome_reference',help="FILENAME genepred format transcriptome")
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('--IDP_isoform_genepred',help='FILENAME IDP output genepred isoform.gpd')
  group1.add_argument('--IDP_output_folder',help="DIRECTORY IDP output directory")
  #parser.add_argument('idp_output_folder',help="FOLDERNAME IDP output")
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--tempdir',default='/tmp',help="FOLDERNAME temproary directory")
  group2.add_argument('--specific_tempdir',help="FOLDERNAME specific temproary directory")
  parser.add_argument('--threads',default=0,type=int,help="INT thread count default to cpu_count")
  parser.add_argument('--max_exon_length',type=int,default=1700,help="INT exons longer than this will not be considered when guaging total gene length")
  parser.add_argument('--min_transcript_expression',type=float,default=0,help="FLOAT dont consider IDP output transcripts with less than or equal to this expression")
  parser.add_argument('--fpkm_cutoff',type=float,default=1,help="FLOAT this is the minimum standard for on RPKM or FPKM in cufflinks or IDP genes")
  parser.add_argument('--cuff',help="FOLDERNAME cufflinks folder if you've already done it for this transcriptome reference")
  parser.add_argument('--filter_by_IDP_expression',action='store_true',help="require that IDP rpkm also meet the cutoff requirements, by default presence in the IDP genepred is sufficient")
  parser.add_argument('-o','--output',help="FILENAME direct output to a file")
  parser.add_argument('--report',help="FILENAME direct a more raw output of data")
  args = parser.parse_args()
  
  if args.IDP_output_folder:
    args.IDP_output_folder = args.IDP_output_folder.rstrip('/')
  if args.threads == 0:
    args.threads = multiprocessing.cpu_count()
  rnum = randint(1,10000000)
  tdir = ''
  if args.specific_tempdir:
    tdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(tdir):
      os.makedirs(tdir)
  else:
    tdir = args.tempdir.rstrip('/')+'/weirathe.'+str(rnum)
    if not os.path.exists(tdir):
      os.makedirs(tdir)
    #print tdir
  sys.stderr.write("working in "+tdir+"\n")

  #establish output
  of = sys.stdout
  if args.output: of = open(args.output,'w')

  # get our short read expression data
  gene_results = get_cufflinks_expression(args,tdir)

  [detection_results,prediction_results] = get_detection_and_prediction(args)

  # get our gene lengths
  lengths = get_lengths(args,tdir)

  # Remove genes from the analysis when they lack any exons meeting our criteria for inclusion.
  # This should only happen for a few genes with exons exceeding the length of the 
  # threshold in the settings.  IF the gene has a transcript in the allowed parameters
  # it will still be included
  glist = gene_results.keys()
  for gene in glist:  
    if gene not in lengths: del gene_results[gene]
  glist = detection_results.keys()
  for gene in glist:  
    if gene not in lengths: del detection_results[gene]
  glist = prediction_results.keys()
  for gene in glist:  
    if gene not in lengths: del prediction_results[gene]
  
  # establish bins
  bins = []
  binsize = 500
  binstop = 5500
  for i in range(0,binstop,binsize):
    e = [i,i+binsize]
    measurements = {}
    measurements['reference'] = set()
    measurements['detected'] = 0
    measurements['predicted'] = 0
    bins.append([e,measurements])
  bins[-1][0][1] = 1000000000
  gene_size_category = {}
  for bin_things in bins:
    genes = get_genes(lengths,gene_results,gene_results,bin_things[0],args.fpkm_cutoff,True)
    for gene in genes: gene_size_category[gene]=str(bin_things[0])
    detected_genes=get_genes(lengths,gene_results,detection_results,bin_things[0],args.fpkm_cutoff,args.filter_by_IDP_expression)
    prediction_genes=get_genes(lengths,gene_results,prediction_results,bin_things[0],args.fpkm_cutoff,args.filter_by_IDP_expression)
    of.write(str(bin_things[0]) + "\t" + str(len(genes)) + "\t" + str(len(detected_genes)) + "\t" + str(len(prediction_genes))+"\n")
    #print str(bin_things[0]) + "\t" + str(genes) + "\t" + str(detected_genes) + "\t" + str(prediction_genes)

  

  if args.report:
    output_report(args,gene_results,gene_size_category,lengths,detection_results,prediction_results)

  if not args.specific_tempdir:
    rmtree(tdir)

def output_report(args,gene_results,gene_size_category,lengths,detection_results,prediction_results):
  of = open(args.report,'w')
  of.write("#Gene\tLength\tLength_Category\tCufflinks_FPKM\tIDP_RPKM\tIDP_Observed\tIDP_Detected\tIDP_Predicted\n")
  for gene in sorted(gene_results.keys()):
    if gene_results[gene] < args.fpkm_cutoff: continue
    scat = ''
    if gene in gene_size_category: scat = gene_size_category[gene]
    express = ''
    dres = 'N'
    if gene in detection_results: 
      dres = 'Y'
      express = str(detection_results[gene])
    pres = 'N'
    if gene in prediction_results:
      pres = 'Y'
      express = str(prediction_results[gene])
    found = 'N'
    if pres == 'Y' or dres == 'Y': found = 'Y'
    of.write(gene +"\t"+str(lengths[gene])+"\t"+scat+"\t" + str(gene_results[gene]) + "\t"+ express+"\t" + found + "\t" + dres + "\t" + pres+"\n")
  of.close()
  
def get_detection_and_prediction(args):
  #get the gene name for every transcript
  conv = {}
  fname = ''
  if args.IDP_output_folder: fname = args.IDP_output_folder+'/isoform.gpd'
  else: fname = args.IDP_isoform_genepred
  with open(fname) as inf:
    for line in inf:
      if re.match('^#',line): continue
      f = line.rstrip().split("\t")
      conv[f[1]] = f[0]
  if args.IDP_output_folder:
    return use_idp_expression(conv,args)
  else:
    return do_not_use_idp_expression(conv,args)

def do_not_use_idp_expression(conv,args):
  # get the expression level of all genes
  predict_gene_expression = {}
  detect_gene_expression = {}
  for transcript in conv:
    if re.search(':\d+-\d+',transcript):
      predict_gene_expression[conv[transcript]] = -1
    else:
      detect_gene_expression[conv[transcript]] = -1
  # remove detection results from prediction
  for gene in detect_gene_expression:
    if gene in predict_gene_expression: del predict_gene_expression[gene]
  return [detect_gene_expression,predict_gene_expression]

def use_idp_expression(conv,args):
  # get the expression level of all genes
  predict_gene_expression = {}
  detect_gene_expression = {}
  with open(args.IDP_output_folder+'/isoform.exp') as inf:
    for line in inf:
      if re.match('^#',line): continue
      f = line.rstrip().split("\t")
      transcript = f[0]
      tx_exp = float(f[1])
      gene_exp = float(f[2])
      if tx_exp <= args.min_transcript_expression and args.min_transcript_expression:
        continue
      if re.search(':\d+-\d+',transcript):
        predict_gene_expression[conv[transcript]] = gene_exp
      else:
        detect_gene_expression[conv[transcript]] = gene_exp
  # remove detection results from prediction
  for gene in detect_gene_expression:
    if gene in predict_gene_expression: del predict_gene_expression[gene]
  return [detect_gene_expression,predict_gene_expression]

# There are two filtering steps here.  Only output results in filter_results passing fpkm
# and then if it passes that if gene_results passes fpkm, then you can do it
# and it must be in lengths
def get_genes(lengths,filter_results,gene_results,bin,fpkm,use_gene_results_fpkm_filter):
  res = set()
  for gene in lengths:
    if lengths[gene] > bin[0] and lengths[gene] <= bin[1]:
      if gene in gene_results:
        # be true if we aren't using this filter (default case)
        if gene_results[gene] >= fpkm or not use_gene_results_fpkm_filter:
          #filter on filter_results
          if gene in filter_results:
            if filter_results[gene] >= fpkm:
              res.add(gene)
  return res

def get_lengths(args,tdir):
  # get our gene lengths
  lengths = {}
  of = open(tdir+'/ref.bed','w')
  cmd = 'bedtools sort -i - | bedtools merge -i - -c 4 -o collapse'
  p = subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=of)
  with open(args.transcriptome_reference) as inf:
    for line in inf:
      if re.match('^#',line): continue
      e = genepred_line_to_entry(line)
      for i in range(0,len(e['exonStarts'])):
        #dont' consider exons that are too long.
        if e['exonEnds'][i]-e['exonStarts'][i] > args.max_exon_length and args.max_exon_length > 0:
          continue
        p.stdin.write(e['chrom']+"\t"+str(e['exonStarts'][i])+"\t"+str(e['exonEnds'][i])+"\t"+e['gene_name']+"\n")
  p.communicate()
  of.close()
  with open(tdir+'/ref.bed') as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      elen = int(f[2])-int(f[1])
      genes = f[3].split(',')
      for gene in genes:
        if gene not in lengths:
          lengths[gene] = 0
        lengths[gene] += elen
  return lengths
      
def get_cufflinks_expression(args,tdir):
  if args.cuff:
    copytree(args.cuff,tdir+'/cuff')
  else:
    # Work on short read alignment by cufflinks
    # 1. make our genepred transcriptome into a gtf
    cmd1 = 'gpd_to_gtf_exons.py '+args.transcriptome_reference+' > '+tdir+'/ref.gtf'
    subprocess.call(cmd1,shell=True)
    # 2. make run our cufflinks
    cmd2 = 'cufflinks -p '+str(args.threads)+' -G  '+tdir+'/ref.gtf '+args.short_read_alignment+' -o '+tdir+'/cuff'
    subprocess.call(cmd2,shell=True)
  
  gene_results = {}
  #Parse cufflinks output
  with open(tdir+'/cuff/genes.fpkm_tracking') as inf:
    header = inf.readline()
    for line in inf:
      f = line.rstrip().split("\t")
      gene_results[f[0]] = float(f[9])
  
  return gene_results

if __name__=="__main__":
  main()
