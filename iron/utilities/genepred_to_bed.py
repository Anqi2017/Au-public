#!/usr/bin/python

import sys, re, os, inspect

#pre: genepred filename, an optional integer for smoothing size
#     The smoothing is done to combine adjacent "exons" within that size definition to get rid of small deletions or indels
#post: a bed file where each line has chromosome, index 1 start, and index 1 end coordinates and the gene_name from each psl line is listed

#bring in the folder to the path for our modules
#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)
import genepred_basics


def main():
  if len(sys.argv) < 2:
    print sys.argv[0] + ' <genepred filename> <smoothing size (optional)>'
    sys.exit()
  genepredfilename = sys.argv[1]
  smooth_size = 0
  if len(sys.argv) == 3:
    smooth_size = int(sys.argv[2])
  with open(genepredfilename) as infile:
    for line in infile:
      if re.match('^#',line):
        continue
      genepred_entry = genepred_basics.genepred_line_to_dictionary(line)
      if smooth_size > 0:
        genepred_entry = genepred_basics.smooth_gaps(genepred_entry,smooth_size)
      for i in range(0,len(genepred_entry['exonStarts'])):
        print genepred_entry['chrom'] + "\t" + str(genepred_entry['exonStarts'][i]+1) + "\t" + str(genepred_entry['exonEnds'][i]) + "\t" + genepred_entry['gene_name'] + "\t" + genepred_entry['name']
main()
