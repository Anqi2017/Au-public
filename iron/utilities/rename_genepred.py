#!/usr/bin/python
import sys, argparse
from GenePredBasics import GenePredEntry
from RangeBasics import Loci, Locus

def main():
  parser = argparse.ArgumentParser(description="Rename gene and transcript elements of GenePred file that are redundant.  Please specify an output if you would like report files generated for the filters.")
  parser.add_argument('input',help="GENEPREDFILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUT FILE default is STDOUT, but you need to specify an output file to get report files generated")
  parser.add_argument('--minimum_locus_distance',type=int,default=500000,help="Genes with the same name will be renamed if this far apart")
  parser.add_argument('--keep_positional_duplicates',action='store_true',help="By default we remove one of the duplicate entries")
  parser.add_argument('--keep_transcript_names',action='store_true',help="By default we rename duplicated transcript names")
  parser.add_argument('--keep_gene_names',action='store_true',help="By default we rename genes located at different loci.")
  args = parser.parse_args()
  inf = sys.stdin
  if args.input != '-': inf = open(args.input)
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  txdef = {}
  gfams = {}
  for line in inf:
    if line[0] == '#': continue
    g = GenePredEntry(line)
    loc = g.value('chrom') + ':' +','.join([str(x) for x in g.value('exonStarts')]) + '-' + ','.join([str(x) for x in g.value('exonEnds')])+'/'+g.value('strand')
    if loc not in txdef:
      txdef[loc] = []
    txdef[loc].append(g)
    if g.value('gene_name') not in gfams: gfams[g.value('gene_name')] = []
    gfams[g.value('gene_name')].append(g.value('name'))
  # now we have cataloged all transcripts by unique locations
  omissions = []
  keepers = []
  for loc in sorted(txdef.keys()):
    if args.keep_positional_duplicates: # We don't want to ommit anything here
      for g in txdef[loc]: keepers.append(g)
      continue #basically skipping this part by populating keepers
    num = len(txdef[loc])
    if num > 1:
      sys.stderr.write("Found "+str(num)+" entries at location\n")
      sys.stderr.write(loc +"\n")
      sys.stderr.write("They are:\n")
      largest = 0
      keepgene = None
      keepindex = -1
      i = 0
      for e in txdef[loc]:
        famsize = len(gfams[e.value('gene_name')])
        sys.stderr.write("     "+e.value('gene_name')+"\t"+e.value('name')+"\t"+str(famsize)+"\n")
        if famsize > largest:
          keepgene = e
          largest = famsize
          keepindex = i
        i+=1
      for j in range(0,len(txdef[loc])):  
        if j != keepindex: omissions.append(txdef[loc][j])
        else: keepers.append(txdef[loc][j])
      sys.stderr.write("     Biggest gene family is "+keepgene.value('gene_name')+" with "+str(largest)+" transcripts\n")
      sys.stderr.write("     so keep that one.\n")
    else:
      keepers.append(txdef[loc][0])
  sys.stderr.write("Omitting "+str(len(omissions))+" entries for redundant positions\n")
  if args.output and not args.keep_positional_duplicates:
    of1 = open(args.output+'.positional_duplicate_omissions','w')
    for g in omissions:
      of1.write(g.get_line()+"\n")
    of1.close()
  # Now the keepers contains transcripts with unique locations
  # Lets provide unique names to remaining transcripts
  tnames = {}
  renametx = {}
  for g in keepers:
    tx = g.value('name')
    if tx not in tnames: tnames[tx] = []
    tnames[tx].append(g)
  for name in tnames:
    if args.keep_transcript_names: continue # We don't want to rename them
    nsize = len(tnames[name])
    if nsize > 1:
      sys.stderr.write("Name: "+name+" has a family of size "+str(nsize)+"\n")
      for i in range(0,len(tnames[name])):
        newname = name+'['+str(i+1)+'/'+str(nsize)+']'
        renametx[newname] = name
        tnames[name][i].entry['name'] = newname
  sys.stderr.write("Renamed: "+str(len(renametx))+" transcripts\n")
  if args.output and not args.keep_transcript_names:
    of1 = open(args.output+'.renamed_transcripts','w')
    for name in sorted(renametx.keys()):
      of1.write(name+"\t"+renametx[name]+"\n")
    of1.close()
  #now we need to arrange into gene families
  gnames = {}
  for name in tnames:
    for g in tnames[name]:
      gene = g.value('gene_name')
      if gene not in gnames:  gnames[gene] = []
      gnames[gene].append(g)
  renamegene = {}
  finished = []
  for gene in gnames:
    if args.keep_gene_names:
      for g in gnames[gene]: finished.append(g)
      continue # We don't want to rename genes
    if len(gnames[gene])==1:
      finished.append(gnames[gene][0])
      continue
    # Now we need to make sure these genes are really on the same locus.
    loci = Loci()
    loci.set_minimum_distance(args.minimum_locus_distance)
    for g in gnames[gene]:
      r = g.locus_range.copy()
      r.set_payload(g)
      loc = Locus()
      loc.add_member(r)
      loci.add_locus(loc)
    loci.update_loci()
    lcount = len(loci.loci)
    if lcount == 1:
      for g in gnames[gene]: finished.append(g)
      continue
    # need to rename some genes
    for i in range(0,lcount):
      newname = gene+'['+str(i+1)+'/'+str(lcount)+']'
      rstr = loci.loci[i].range.get_range_string()
      renamegene[newname] = gene
      sys.stderr.write(newname+"\t"+rstr+"\n")
      for m in loci.loci[i].members:
        m.get_payload().entry['gene_name'] = newname
        finished.append(m.get_payload())
  sys.stderr.write("Renamed: "+str(len(renamegene))+" genes\n")
  if args.output and not args.keep_transcript_names:
    of1 = open(args.output+'.renamed_genes','w')
    for name in sorted(renamegene.keys()):
      of1.write(name+"\t"+renamegene[name]+"\n")
    of1.close()
  #Now lets resort by genes
  bygene = {}
  for g in finished:
    gene = g.value('gene_name')
    if gene not in bygene: bygene[gene] = []
    bygene[gene].append(g)
  for gene in sorted(bygene.keys()):
    for g in bygene[gene]:
      of.write(g.get_line()+"\n")
  of.close()
  inf.close()

if __name__=="__main__":
  main()
