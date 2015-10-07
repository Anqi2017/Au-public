#!/usr/bin/python
import argparse, re, sys, multiprocessing, json
import PSLBasics, GenePredBasics
from FileBasics import GenericFileReader
from SequenceBasics import read_fasta_into_hash

local_support = {}

def main():
  parser = argparse.ArgumentParser(description='Use reference junctions when they are close',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--min_intron_size',type=int,default=68,help="INT min intron size")
  parser.add_argument('--min_local_support',type=int,default=0,help="INT min number of junctions within search_size of a junction in order to count it")
  parser.add_argument('--search_size',type=int,default=10,help="INT search space for reference")
  parser.add_argument('--output_fake_psl',help="FASTAFILE reference genome to make a fake PSL output")
  parser.add_argument('psl',help="PSLFILENAME or '-' for STDIN")
  parser.add_argument('reference_genepred',help="FASTAFILENAME for reference genepred")
  args = parser.parse_args()

  cpus = multiprocessing.cpu_count()

  genome = {}
  if args.output_fake_psl:
    genome = read_fasta_into_hash(args.output_fake_psl)

  #read in the reference genepred first
  gpf = GenePredBasics.GenePredFile(args.reference_genepred)
  #lets sort entries by chromosome
  ref = {}
  for e in [x.entry for x in gpf.entries]:
    if len(e['exonStarts']) <= 1: continue
    if e['chrom'] not in ref:
      ref[e['chrom']] = {}
    for i in range(1,len(e['exonStarts'])):
      if e['exonEnds'][i-1] not in ref[e['chrom']]:
        ref[e['chrom']][e['exonEnds'][i-1]] = {}
      if e['exonStarts'][i]+1 not in ref[e['chrom']][e['exonEnds'][i-1]]:
        ref[e['chrom']][e['exonEnds'][i-1]][e['exonStarts'][i]+1] = e['strand']
  #Stored all junctions as 1-base

  read_info = {}
  pf = GenericFileReader(args.psl)
  fcount_total = 0
  while True:
    line = pf.readline()
    if not line: break
    if re.match('^#',line): continue
    line = line.rstrip()
    pe = PSLBasics.line_to_entry(line)
    if len(pe['tStarts']) != len(pe['blockSizes']) or len(pe['qStarts']) != len(pe['blockSizes']):
      sys.stderr.write("WARNING invalid psl\n")
      continue
    genepred_line = PSLBasics.convert_entry_to_genepred_line(pe)
    ge = GenePredBasics.smooth_gaps(GenePredBasics.line_to_entry(genepred_line),args.min_intron_size)
    refjuns = {}
    if pe['tName'] in ref: refjuns = ref[pe['tName']]
    new_ge = nudge(pe,ge,refjuns,args)
    if args.output_fake_psl:
      new_psl_line = GenePredBasics.entry_to_fake_psl_line(new_ge,genome)
      print new_psl_line
    else:
      print GenePredBasics.entry_to_line(new_ge)

def nudge(psl_entry,gpd_entry,refjun,args):
  junctions = []
  fcount = 0
  if len(gpd_entry['exonStarts']) == 1:
    #print "no intron 1"
    return gpd_entry
  bounds = []
  for i in range(1,len(gpd_entry['exonStarts'])):
    junc_start = gpd_entry['exonEnds'][i-1]
    junc_finish = gpd_entry['exonStarts'][i]+1
    bounds.append([junc_start, junc_finish,i-1])
  if len(bounds) < 1:
    #print "no intron 2"
    return gpd_entry
  bestbounds = []
  for bound in bounds:
    best_distance = [10000000,10000000]
    best_result = None
    for z1 in range(bound[0]-args.search_size,bound[0]+args.search_size+1):
      d1 = abs(z1-bound[0])
      if z1 in refjun:
        for z2 in range(bound[1]-args.search_size,bound[1]+args.search_size+args.search_size+1):
          d2 = abs(z2-bound[1])
          if z2 in refjun[z1]:
            refstrand = refjun[z1][z2]
            if d1+d2 < best_distance[0]+best_distance[1]:
              best_distance = [d1,d2]
              best_result = [z1,z2,refstrand,bound[2]]+best_distance
    if best_result:
      bestbounds.append(best_result)
  if len(bestbounds) < 1: 
    #nothing fixable
    #sys.stderr.write("nothing fixable\n")
    return gpd_entry
  #Now we have a list of nudgable bounds
  #Lets pick a strand
  plus_score = 0
  minus_score = 0
  #print '----'
  #print bestbounds
  for bound in bestbounds:
    if bound[2] == '+':
      plus_score += 1/(float(abs(bound[4]))+float(abs(bound[5]))+1)
    else:
      minus_score += 1/(float(abs(bound[4]))+float(abs(bound[5]))+1)
  use_strand = '+'
  #print [plus_score,minus_score]
  if plus_score < minus_score: use_strand = '-'
  #print use_strand
  choice_bounds = []
  for bound in bestbounds:
    if bound[2] == use_strand:  choice_bounds.append(bound)
  #print '---'
  #print GenePredBasics.entry_to_line(gpd_entry)
  #print bestbounds
  #print choice_bounds
  if len(choice_bounds) < 1: 
    print "ERROR  should have choices"
    sys.exit()
  replacements = {}
  for bound in choice_bounds:  replacements[bound[3]] = [bound[0],bound[1]]
  junctions = []
  #print "fixed "+str(len(replacements.keys()))
  for i in range(0,len(bounds)):
    val = bounds[i]
    if i in replacements:
      #sys.stderr.write("use replacement\n")
      val = replacements[i]
      fcount += 1
    junctions.append([val[0],val[1]])
  #print junctions
  #sys.stderr.write("replace\n")
  #print junctions
  new_gpd_line  = gpd_entry['gene_name'] + "\t"
  new_gpd_line += gpd_entry['name'] + "\t"
  new_gpd_line += gpd_entry['chrom'] + "\t"
  new_gpd_line += gpd_entry['strand'] + "\t"
  new_gpd_line += str(gpd_entry['txStart']) + "\t"
  new_gpd_line += str(gpd_entry['txEnd']) + "\t"
  new_gpd_line += str(gpd_entry['cdsStart']) + "\t"
  new_gpd_line += str(gpd_entry['cdsEnd']) + "\t"
  new_gpd_line += str(len(junctions)+1) + "\t"
  exon_starts = [gpd_entry['txStart']]
  exon_ends = [] #gpd_entry['txEnd']]
  for junc in junctions:
    exon_starts.append(junc[1]-1)
    exon_ends.append(junc[0])
  exon_ends.append(gpd_entry['txEnd'])
  new_gpd_line += ','.join([str(x) for x in exon_starts])+','+"\t"
  new_gpd_line += ','.join([str(x) for x in exon_ends])+','+"\t"
  #print new_gpd_line
  new_gpd_entry = GenePredBasics.line_to_entry(new_gpd_line)
  #print "got junctions"
  #print new_gpd_line
  #print '.........'
  return new_gpd_entry
  
if __name__=="__main__":
  main()
