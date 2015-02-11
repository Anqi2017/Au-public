#!/usr/bin/python

import sys, os, inspect, re
from copy import deepcopy
from shutil import rmtree, copyfile

#Pre: <long reads fasta> <reference genome fasta> <refFlat genepred> <out base>
#     Requires gmap and bedtools be installed in the path
#Post: 
#Modifies: Requires a nonrepetative genepred file.  only one entry per transcript
#          Only one entry per set of coordinates.  If this doesnt exist, it will
#          make it and use it, but it is greedy and will only do the first
#          entry of each repeat name or coordinate
#          by name i mean the second field in the genepred file which should
#          be the transcript name

#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import aligner_basics, genepred_basics, FileBasics, psl_basics


def main():
  if len(sys.argv) != 5:
    print sys.argv[0]+' <long reads fasta> <reference genome> <transcriptome genepred>'
    return
  longreadfname = sys.argv[1]
  genomefname = sys.argv[2]
  usergenepredfname = sys.argv[3]
  outbase = sys.argv[4]

  #get read count
  readcount = 0
  with open(longreadfname) as f:
    for line in f:
      if re.match('^>',line): readcount+=1

  tdir = FileBasics.make_tempdir2('weirathe','annlong')

  # 1. Make sure the transcriptome is uniquely named uniquely mapped entries
  genepredfname = tdir+'/txn.gpd'
  make_unique_genepred(usergenepredfname,genepredfname)
  print 'made unique genepred file'

  # 2.  Make a transcriptome to align to.
  transcriptomefasta = tdir+'/txn.fa'
  genepred_basics.write_genepred_to_fasta_directionless(genepredfname,genomefname,transcriptomefasta)
  print 'made transcriptome fasta'

  # 3.  Make a bed file of junction locations in that transcriptome.
  junctionbedfname = tdir+'/junction.bed'
  junction_counts = make_junction_bed_file(genepredfname,junctionbedfname)
  print 'made junction bed file'

  # 4.  Build a gmap index of the transcriptome
  transcriptomeindex = tdir+'/gmap_txn'
  aligner_basics.build_gmap_index(transcriptomefasta,transcriptomeindex)
  print 'made gmap index of transcriptome'

  # 5.  Align the long reads to the transcriptome with gmap
  alignmentfname = tdir+'/reads.psl'
  aligner_basics.gmap_all(longreadfname,transcriptomeindex,alignmentfname)
  print 'made gmap alignment of reads to transcriptome'

  # 6.  Generate get the genepred of the long reads on the transcriptome coordinates.
  #     Smooth that genepred by a smoothing factor
  #     And make a bed file of the best alignment.  see function for specifications
  bestalignmentbedfname = tdir+'/reads.bed'
  make_best_continuous_alignment_bed(alignmentfname,bestalignmentbedfname)
  print 'made best continuous alignment bed file'

  # 10.  Print per-gene count info
  genenames = genepred_basics.get_transcript_to_gene_name_dictionary(genepredfname)
  print 'got gene name conversions'  

  # 7.  Make a report of all prefilter alignments  
  bestprefilter = tdir+'/prefilter.txt'
  prefilter_alignments = make_best_alignment_summary(bestalignmentbedfname,junctionbedfname,junction_counts,bestprefilter)
  print 'made best alignment prefilter summary'

  report_file = tdir +'/report.txt'
  orep = open(report_file,'w')
  orep.write('Basename:'+"\t"+outbase+"\n")
  orep.write('Temp directory:'+"\t"+tdir+"\n")
  orep.write('Long Read Count:'+"\t"+str(readcount)+"\n")
  # 8.  Filter the full length alignments 
  full_length_alignments = filter_alignments(prefilter_alignments,'full')
  full_length_alignment_file = tdir+'/full_length_alignment.txt'
  [full_length_read_count, full_length_transcript_count] = write_alignments(full_length_alignments,full_length_alignment_file,genenames)
  orep.write('Read count - full length reads mapped:'+"\t" +str(len(full_length_alignments))+"\n")
  orep.write('Transcript count - full length reads mapped:'+"\t" +str(full_length_transcript_count)+"\n")
  unambiguous_full_length_alignment_file = tdir+'/unambiguous_full_length_alignment.txt'
  unambiguous_full_length_alignments = filter_unambiguous_alignments(full_length_alignments)
  [unambiguous_full_length_read_count, unambiguous_full_length_transcript_count] = write_alignments(unambiguous_full_length_alignments,unambiguous_full_length_alignment_file,genenames)
  orep.write('Read count - full length reads mapped with unambiguous matches:'+"\t"+str(len(unambiguous_full_length_alignments))+"\n")
  orep.write('Transcript count - full length reads mapped with unambiguous matches:'+"\t"+str(unambiguous_full_length_transcript_count)+"\n")

  # 9.  Filter the full length alignments 
  prepartial_alignments = filter_alignments(prefilter_alignments,'partial')
  prepartial_alignment_file = tdir+'/prepartial_alignment.txt'
  write_alignments(prepartial_alignments,prepartial_alignment_file,genenames)
  partial_alignments = filter_by_priority_alignments(prepartial_alignments)
  partial_alignment_file = tdir+'/partial_alignment.txt'
  [partial_read_count, partial_transcript_count] = write_alignments(partial_alignments,partial_alignment_file,genenames)
  orep.write('Read count - reads mapped with partial hits best junction and length matches:'+"\t" + str(len(partial_alignments))+"\n")
  orep.write('Transcript count - reads mapped with partial hits best junction and length matches:'+"\t" + str(partial_transcript_count)+"\n")
  unambiguous_partial_alignments = filter_unambiguous_alignments(partial_alignments)
  unambiguous_partial_alignment_file = tdir + '/unambiguous_partial_alignments.txt'
  [unambiguous_partial_read_count, unambiguous_partial_transcript_count] = write_alignments(unambiguous_partial_alignments,unambiguous_partial_alignment_file,genenames)
  orep.write('Read count - reads mapped with partial hits unambiguous matches:'+"\t"+str(len(unambiguous_partial_alignments))+"\n")
  orep.write('Transcript count - reads mapped with partial hits unambiguous matches:'+"\t"+str(unambiguous_partial_transcript_count)+"\n")

  partial_gene_counts = get_uniquely_mappable_gene_counts(partial_alignments,genenames)
  partial_gene_counts_file = tdir+'/partial_match_uniquely_mappable_gene_counts.txt'
  write_gene_counts(partial_gene_counts,partial_gene_counts_file)

  full_gene_counts = get_uniquely_mappable_gene_counts(full_length_alignments,genenames)
  full_gene_counts_file = tdir+'/full_length_match_uniquely_mappable_gene_counts.txt'
  write_gene_counts(full_gene_counts,full_gene_counts_file)

  orep.write('Gene count - full length matches uniquely mapped:'+"\t"+str(len(full_gene_counts))+"\n")
  orep.write('Gene count - partial matches uniquely mapped:'+"\t"+str(len(partial_gene_counts))+"\n")
  orep.close()
  copyfile(report_file,outbase+'.Report.txt')
  copyfile(full_gene_counts_file,outbase+'.FullGeneCounts.txt')
  copyfile(partial_gene_counts_file,outbase+'.PartialGeneCounts.txt')
  copyfile(full_length_alignment_file,outbase+'.FullAlignment.txt')
  copyfile(unambiguous_full_length_alignment_file,outbase+'.UnambiguousFullAlignment.txt')
  copyfile(partial_alignment_file,outbase+'.PartialAlignment.txt')
  copyfile(unambiguous_partial_alignment_file,outbase+'.UnambiguousFullAlignment.txt')
  rmtree(tdir)

#Post: <transcript> <gene> <read> <number of transcripts this read maches to>
#      <junctions spanned> <junctions in the transcript> <length of alignment>
#      <length of read> <length of transcript> <number of matches>
def write_alignments(alignments,full_length_alignment_file,conv):
  of = open(full_length_alignment_file,'w')
  rds = set()
  txs = set()
  for read in alignments:
    rds.add(read)
    num = len(alignments[read])
    for tx in alignments[read]:
      txs.add(tx)
      e = alignments[read][tx]
      gene = conv[tx]
      of.write(tx + "\t" + gene + "\t" + read + "\t" + str(num) + "\t" + \
            str(e['junctions_seen_count']) + "\t" + str(e['total_junctions']) + "\t" + \
            str(e['alignment_length']) + "\t" + str(e['read_length']) + "\t" + \
            str(e['transcript_length']) + "\t" + str(e['matches']) + "\n")
  of.close()
  return [len(rds),len(txs)]

def filter_unambiguous_alignments(alignment):
  out = {}
  for read in alignment:
    if len(alignment[read]) != 1: continue
    if read not in out:
      out[read] = {}
    for tx in alignment[read]:
      out[read][tx] = deepcopy(alignment[read][tx])
  return out

def write_gene_counts(gene_counts,outfile):
  of = open(outfile,'w')
  for gene in sorted(gene_counts,key=gene_counts.get,reverse=True):
    of.write(gene + "\t" + str(gene_counts[gene])+"\n")
  of.close()

def get_uniquely_mappable_gene_counts(alignment,conv):
  counts = {}
  for read in alignment:
    genes = set()
    for tx in alignment[read]:
      e = alignment[read][tx]
      genename = conv[tx]
      genes.add(genename)
    if len(genes) == 1:
      gname = genes.pop()
      if gname not in counts:
        counts[gname] = 0
      counts[gname] += 1
  return counts


# Pre: each entry keyed by a read
#        e['transcript'] 
#        e['read_length']
#        e['transcript_length']
#        e['alignment_length']
#        e['total_junctions']
#        e['junctions_seen_count']
# Post: Filters by the best junction agreement first
#       then longest length of alignment second
def filter_by_priority_alignments(p):
  out = {}
  for read in p:
    # sort by junctions take the one with the most junctions without too many.
    # get best junction count
    bestjcount = -1
    for tx in p[read]:
      e = p[read][tx]
      #print tx + "\t" +  str(e['total_junctions']) + "\t" +  str(e['alignment_length']) 
      if e['junctions_seen_count'] > bestjcount and e['total_junctions'] >= e['junctions_seen_count']:
        bestjcount = p[read][tx]['junctions_seen_count']
    bestlen = 0
    for tx in p[read]:
      e = p[read][tx]
      if e['alignment_length'] > bestlen and e['junctions_seen_count'] == bestjcount:
        bestlen = e['alignment_length']
    for tx in p[read]:
      e = p[read][tx]
      if e['alignment_length'] == bestlen and e['junctions_seen_count']:
        if read not in out:
          out[read] = {}
        out[read][tx] = deepcopy(p[read][tx])
  return out
    
   

# Pre: each entry keyed by a read
#        e['transcript'] 
#        e['read_length']
#        e['transcript_length']
#        e['alignment_length']
#        e['total_junctions']
#        e['junctions_seen_count']

def filter_alignments(p,type):
  out = {}
  for read in p:
    for tx in p[read]:
      read_cov_by_alignment = 0.8
      min_length_of_alignment = 100
      if p[read][tx]['transcript_length'] == 0:
        sys.stderr.write("strange transcript length of 0. skipping\n")
        continue
      if p[read][tx]['alignment_length'] == 0:
        sys.stderr.write("strange read length of 0 skipping\n")
        continue
      if type == 'full' and p[read][tx]['junctions_seen_count'] != p[read][tx]['total_junctions']:
        continue
      if min_length_of_alignment > p[read][tx]['alignment_length']:
        continue
      overlap = float(p[read][tx]['alignment_length'])/float(p[read][tx]['read_length'])
      if overlap < read_cov_by_alignment:
        continue
      if read not in out:
        out[read] = {}
      out[read][tx] = deepcopy(p[read][tx])
  return out

def make_best_alignment_summary(bestalignmentbedfname,junctionbedfname,junction_counts,bestalignmentsummary):
  of = open(bestalignmentsummary,'w')
  cmd = 'bedtools intersect -wao -a '+bestalignmentbedfname+' -b '+junctionbedfname
  #sys.stderr.write(cmd+"\n")
  reads = {}
  with os.popen(cmd) as pipe:
    for line in pipe:
      f = line.rstrip().split("\t")
      read = f[3]
      transcript = f[0]
      if read not in reads:
        reads[read] = {}
      if transcript not in reads[read]:
        e = {}
        e['read_length'] = int(f[4])
        e['transcript_length'] = int(f[5])
        e['alignment_length'] = int(f[6])
        e['matches'] = int(f[7])
        e['total_junctions'] = junction_counts[transcript]
        e['junctions_seen_count'] = 0
        e['junctions_seen'] = set()
        reads[read][transcript] = e
      junction_overlap = 0      
      if f[12] != '.': 
        junction_overlap = int(f[12])
      if junction_overlap == 2: # left and right base pairs align back
        reads[read][transcript]['junctions_seen'].add(f[11])
  for read in reads:
    for transcript in reads[read]:
      reads[read][transcript]['junctions_seen_count'] = len(reads[read][transcript]['junctions_seen'])
      of.write(read + "\t" + transcript + "\t" + \
               str(reads[read][transcript]['read_length']) + "\t" + str(reads[read][transcript]['transcript_length']) + "\t" + \
               str(reads[read][transcript]['alignment_length']) + "\t" + str(reads[read][transcript]['total_junctions']) + "\t" + \
               str(reads[read][transcript]['junctions_seen_count'])+"\n")
  return reads

#Pre: A psl file of long reads aligned to transcriptome
#Post: A bed file where each entry is the longest continuous alignment
#      of the long reads (query) with the transcriptome (target) in the 
#      following format:
#      <transcript (target)> <start> <end> <read (query)> <read length> <transcript length> <length of aligned region>

def make_best_continuous_alignment_bed(alignmentfname,bestalignmentbedfname):
  bestcont = {}
  of = open(bestalignmentbedfname,'w')
  with open(alignmentfname) as f:
    for line in f:
      if re.match('^#',line): continue
      pe = psl_basics.read_psl_entry(line)
      gl = psl_basics.convert_entry_to_genepred_line(pe)
      ge = genepred_basics.genepred_line_to_dictionary(gl)
      ges = genepred_basics.smooth_gaps(ge,10)
      for i in range(0,len(ges['exonStarts'])):
        exonlen = ges['exonEnds'][i]-ges['exonStarts'][i]
        if ges['name'] not in bestcont:
          bestcont[ges['name']] = {}
        if ges['chrom'] not in bestcont[ges['name']]:
          entry = {}
          entry['bestlen'] = 0
          entry['start'] = 0
          entry['end'] = 0
          entry['target_length'] = 0
          entry['query_length'] = 0
          entry['matches'] = 0
          bestcont[ges['name']][ges['chrom']] = entry
        if exonlen > bestcont[ges['name']][ges['chrom']]['bestlen']:
          bestcont[ges['name']][ges['chrom']]['bestlen'] = exonlen
          bestcont[ges['name']][ges['chrom']]['start'] = ges['exonStarts'][i]
          bestcont[ges['name']][ges['chrom']]['end'] = ges['exonEnds'][i]
          bestcont[ges['name']][ges['chrom']]['target_length'] = pe['tSize']
          bestcont[ges['name']][ges['chrom']]['query_length'] = pe['qSize']
          bestcont[ges['name']][ges['chrom']]['matches'] = pe['matches']
  for read in bestcont:
    for tx in bestcont[read]:
      of.write(tx + "\t" + str(bestcont[read][tx]['start']) + "\t" +  str(bestcont[read][tx]['end']) + "\t" + read + "\t" + str(bestcont[read][tx]['query_length']) + "\t" + str(bestcont[read][tx]['target_length']) + "\t"  + str(bestcont[read][tx]['bestlen']) + "\t" + str(bestcont[read][tx]['matches']) + "\n")
  return
  
# Pre: <genepred filename> <junction filename>
# Post: writes bed file of the following format
#       <transcript name> <junction start> <junction end> <junction unique name>
#       returns a dictionary of transcript names with junction counts
# Modifies:

def make_junction_bed_file(gpd_filename,junc_filename):
  jcounts = {}
  of = open(junc_filename,'w')
  with open(gpd_filename) as f:
    for line in f:
      if re.match('^#',line): continue
      d = genepred_basics.genepred_line_to_dictionary(line)
      juncs = []
      totlen = 0
      for i in range(0,len(d['exonStarts'])):
        totlen += d['exonEnds'][i]-d['exonStarts'][i]
        juncs.append(str(totlen-1)+"\t"+str(totlen+1))
      juncs.pop()
      j = 0
      jcounts[d['name']] = len(juncs)
      for junc in juncs:
        j+=1
        of.write(d['name'] + "\t"+junc+"\t"+d['name']+'.'+str(j) + "\n")
  return jcounts

def make_unique_genepred(genepredfname, uniquegenepredfname):
  of = open(uniquegenepredfname,'w')
  names = set()
  coords = set()
  skipcount = 0
  with open(genepredfname) as f:
    for line in f:
      if re.match('^#',line): continue
      d = genepred_basics.genepred_line_to_dictionary(line)
      coord = d['chrom']+':'+','.join(str(x) for x in d['exonStarts'])+'-'+','.join(str(x) for x in d['exonEnds'])
      if d['name'] in names:
        sys.stderr.write("skipping repeat name: "+line.rstrip()+"\n")
        skipcount += 1
        continue
      if coord in coords:
        sys.stderr.write("skipping repeat coord: "+line.rstrip()+"\n")
        skipcount += 1
        continue
      names.add(d['name'])
      coords.add(coord)
      of.write(line.rstrip()+"\n")
  if skipcount > 0:
    sys.stderr.write("skipped a total of "+str(skipcount)+"\n")
  return

main()
