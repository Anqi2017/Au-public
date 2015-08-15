#!/usr/bin/python
import argparse, sys, re, os
import GenePredBasics

def main():
  parser = argparse.ArgumentParser(description="Make a universal genepred and key for comparing IDP results")
  parser.add_argument('--output_directory',default='IDP_output_merge',help='DIRECTORY to write output to.  Will not overwrite existing')
  parser.add_argument('genepred_exp_name_sets',nargs='+',help="three files for each IDP entry 1) a genepred 2) an expression file 3) a sample name.")
  args = parser.parse_args()

  mydir = args.output_directory.rstrip('/')
  if os.path.isdir(mydir):
    sys.stderr.write("ERROR: output directory "+mydir+" already exists\n")
    return
  os.makedirs(mydir)

  set_args = args.genepred_exp_name_sets
  if len(set_args)%3 != 0:
    sys.stderr.write("Data must be in sets of three")
  setnum = 0
  resultnumber = 0
  numbers = {}
  byset = {}
  chromosomes = set()
  established_names = {}
  expression = {}
  sample_names = set()
  while len(set_args) > 0:
    setnum += 1
    gpd = set_args.pop(0)
    exp = set_args.pop(0)
    sample_name = set_args.pop(0)
    sample_names.add(sample_name)
    sys.stderr.write("Set: "+str(setnum)+"\n")
    sys.stderr.write("  GenePred: "+gpd+"\n")
    sys.stderr.write("  Expression: "+exp+"\n")
    sys.stderr.write("  Sample: "+sample_name+"\n")
    
    with open(gpd) as inf:
      for line in inf:
        if re.match('^#',line): continue
        e = GenePredBasics.GenePredEntry()
        e.line_to_entry(line)
        chromosomes.add(e.entry['chrom'])
        junctions =  e.junctions
        resultnumber += 1
        junstring = ";".join(junctions)
        if junstring not in byset:
          byset[junstring] = set()
        byset[junstring].add(resultnumber)
        numbers[resultnumber] = [sample_name,e.entry['name'],e]
    with open(exp) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        if sample_name not in expression:
          expression[sample_name] = {}
        expression[sample_name][f[0]] = [float(f[1]),float(f[2])] #transcript and gene exression

  #bysample = {}
  gene_records = {}
  for junc in byset:
    lowest = False
    highest = False
    realnames = set()
    realgenenames = set()
    chromnames = set()
    chromgenenames = set()
    arbitrary_gpd = False
    sgpds = {}
    for i in byset[junc]:
      [sample, name,gpd] = numbers[i]
      gene_name = gpd.entry['gene_name']
      arbitrary_gpd = gpd
      sgpds[sample] = gpd
      # Figure out if its a reference transcript name or an IDP manufactured name
      m = re.match('^([^:]+):\d+-\d+',name)
      if not m:
        realnames.add(name)
      else:
        chromnames.add(m.group(1))

      # Figure out if its a reference gene name or an IPD manufacture gene anme
      m = re.match('^([^:]+):\d+-\d+',gene_name)
      if not m:
        realgenenames.add(gene_name)
      else:
        chromgenenames.add(m.group(1))

      if not lowest or gpd.entry['txStart'] < lowest:
        lowest = gpd.entry['txStart']
      if not highest or gpd.entry['txEnd'] > highest:
        highest = gpd.entry['txEnd']
      #if sample not in bysample:
      #  bysample[sample] = {}
      #if name not in bysample[sample]:
      #  bysample[sample][name] = i
    usename = False
    basename = False
    if len(realnames) > 0:
      usename = next(iter(realnames))
      if len(realnames) > 1: 
        sys.stderr.write("WARNING: multiple transcript names as with the same junctions.\n"+str(realnames)+"\nUsing: "+str(usename)+"\n")
      if usename in established_names:
        sys.stderr.write("WARNING: reference transcript name "+usename+" refers to different transcripts with different junction compositions.  Renaming the second instance to a unique name.")
        established_names[usename] += 1
        usename = usename +'.'+str(established_names[usename])
      else:
        established_names[usename] = 0
    else:
      usechrom = next(iter(chromnames))
      if len(chromnames) > 1:
        sys.stderr.write("ERROR: multiple chromosome names are not supported in a single transcript yet.\n"+str(chromnames)+"\n")
        sys.exit()
      basename = usechrom + ":" + str(lowest) + '-' + str(highest)
      if basename not in established_names:
        established_names[basename] = 0
      established_names[basename]+=1
      usename = basename + '.'+str(established_names[basename])
    # See if we have a real gene name for base name
    if len(realgenenames) > 0:
      basename = next(iter(realgenenames))
    #print basename + "\t" + usename 
    if basename not in gene_records:
      gene_records[basename] = {}
    gene_records[basename][usename] = {}
    gene_records[basename][usename]['sample_gpd'] = {}
    gene_records[basename][usename]['sample_exp'] = {}
    gene_records[basename][usename]['gpd'] = GenePredBasics.GenePredEntry()
    # copy the old record
    gene_records[basename][usename]['gpd'].line_to_entry(arbitrary_gpd.get_line())
    if lowest < gene_records[basename][usename]['gpd'].entry['txStart']:
      sys.stderr.write("ADJUSTING NEW GPD TXSTART FOR "+basename+" " + usename+"\n")
      gene_records[basename][usename]['gpd'].entry['txStart'] = lowest
      gene_records[basename][usename]['gpd'].entry['cdsStart'] = lowest
      gene_records[basename][usename]['gpd'].entry['exonStarts'][0] = lowest

    if highest > gene_records[basename][usename]['gpd'].entry['txEnd']:
      sys.stderr.write("ADJUSTING NEW GPD TXEND FOR "+basename+" " + usename+"\n")
      gene_records[basename][usename]['gpd'].entry['txEnd'] = highest
      gene_records[basename][usename]['gpd'].entry['cdsEnd'] = highest
      gene_records[basename][usename]['gpd'].entry['exonEnds'][len(gene_records[basename][usename]['gpd'].entry['exonEnds'])-1] = lowest
    # Now add the original sample information     
    for sample in sgpds:
      gene_records[basename][usename]['sample_gpd'][sample]=sgpds[sample]
      gene_records[basename][usename]['sample_exp'][sample]=expression[sample][sgpds[sample].entry['name']][0]

  #Now all necessary data should be in gene_records
  sample_list = sorted(list(sample_names))
  ofgene = open(mydir+'/gene.exp','w')
  ofgene.write("gene")
  for sample in sample_list:
    ofgene.write("\t"+sample)
  ofgene.write("\n")
  geneexp = {}
  for gene in gene_records:
    total = {}
    for sample in sample_list:  total[sample] = 0
    geneexp[gene] = {}
    for transcript in gene_records[gene]:
      for sample in gene_records[gene][transcript]['sample_exp']:
        total[sample] += gene_records[gene][transcript]['sample_exp'][sample]
    ofgene.write(gene)
    for sample in sample_list:
      geneexp[gene][sample]  = total[sample]
      ofgene.write("\t"+str(total[sample]))
    ofgene.write("\n")
  ofgene.close()

  #Now we can do all the transcript writing
  ofgeneiso = open(mydir+'/gene_isoform.exp','w')
  ofgeneiso.write("gene\tisoform")
  for sample in sample_list:
    ofgeneiso.write("\t"+sample+".gene"+"\t"+sample+".isoform")
  ofgeneiso.write("\n")
  ofiso = open(mydir+'/isoform.exp','w')
  ofiso.write("isoform")
  for sample in sample_list:
    ofiso.write("\t"+sample)
  ofiso.write("\n")
  for gene in gene_records:
    for transcript in gene_records[gene]:
      ofiso.write(transcript)
      ofgeneiso.write(gene+"\t"+transcript)
      for sample in sample_list:
        if sample in gene_records[gene][transcript]['sample_exp']:
          ofgeneiso.write("\t"+str(geneexp[gene][sample])+"\t"+str(gene_records[gene][transcript]['sample_exp'][sample]))
          ofiso.write("\t"+str(gene_records[gene][transcript]['sample_exp'][sample]))
        else:
          ofiso.write("\t0")
          if sample in geneexp[gene]:
            ofgeneiso.write("\t"+str(geneexp[gene][sample])+"\t0")
          else:
            ofgeneiso.write("0\t0")
      ofiso.write("\n")
      ofgeneiso.write("\n")
  ofiso.close()

  #Maybe we can finish it all off by writing the new genepred
  ofgpd = open(mydir+'/isoform.gpd','w')
  for gene in gene_records:
    for transcript in gene_records[gene]:
      ofgpd.write(gene_records[gene][transcript]['gpd'].get_line()+"\n")
  ofgpd.close()

if __name__=="__main__":
  main()
