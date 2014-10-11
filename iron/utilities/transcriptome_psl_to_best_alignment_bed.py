#!/usr/bin/python
import sys, os, inspect
#      This can be used as part of a pipeline to call which transcript
#      and long read is aligned to.  This takes the alignment of the long
#      read to the transcriptome, and it outputs a sit of possible alignments
#      that are in bed format to be more easily filtered and scanned for 
#      desired junction/coverage requirements
#
#pre: a psl file from reads aligned to a directionless transcriptome fasta
#     directionless means that transcripts all go in coordinate order and
#     sequences were not reverse complimented for negative strands
#     see genepred_basics write directionless fasta function
#post:  A bed file with the following best continuous alignment reported for each psl entry
#     <target name> <start 0-indexed> <end 1-indexed> <query name> <query length> <target length> <best continuous alignment length>
#      That means if multiple exons are reported for a single psl entry,
#      only the longest exon will be returned.
#modifies fileIO

#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)


import genepred_basics, psl_basics


def main():

  if len(sys.argv) < 2:
    print sys.argv[0] + " <psl filename> <smoothing parameter (default 10)>"
    sys.exit()
  smoothing = 10
  if len(sys.argv) == 3:
    smoothing = int(sys.argv[2])
  psl_filename = sys.argv[1]
  with open(psl_filename) as fh:
    for line in fh:
      line = line.rstrip()
      psl_entry = psl_basics.read_psl_entry(line)
      gpd_line = psl_basics.convert_entry_to_genepred_line(psl_entry)
      gpd_entry = genepred_basics.genepred_line_to_dictionary(gpd_line)
      smoothed_gpd_entry = genepred_basics.smooth_gaps(gpd_entry,smoothing)
      #get longest exon only
      longest_exon = 0
      best_start = 0
      best_end = 0
      for i in range(0,smoothed_gpd_entry['exonCount']):
        exon_length = smoothed_gpd_entry['exonEnds'][i]-smoothed_gpd_entry['exonStarts'][i]
        if exon_length > longest_exon:
          longest_exon = exon_length
          best_end = smoothed_gpd_entry['exonEnds'][i]
          best_start = smoothed_gpd_entry['exonStarts'][i]
      print psl_entry['tName'] + "\t" +  str(best_start) + "\t" + str(best_end) + "\t" + psl_entry['qName'] + "\t" + str(psl_entry['qSize']) + "\t" + str(psl_entry['tSize']) + "\t" + str(longest_exon)

main()
