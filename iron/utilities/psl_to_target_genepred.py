#!/usr/bin/python

import sys, os, inspect

#pre: psl filename, an optional integer for smoothing size
#     The smoothing is done to combine adjacent "exons" within that size definition to get rid of small deletions or indels
#post: a genepred file with target sequence information.  not query

#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)
import psl_basics, genepred_basics


def main():
  if len(sys.argv) < 2:
    print sys.argv[0] + ' <psl filename> <smoothing size (optional)>'
    sys.exit()
  pslfilename = sys.argv[1]
  smooth_size = 0
  if len(sys.argv) == 3:
    smooth_size = int(sys.argv[2])
  with open(pslfilename) as infile:
    for line in infile:
      psl_entry = psl_basics.read_psl_entry(line)
      genepred_line = psl_basics.convert_entry_to_genepred_line(psl_entry)
      print genepred_line
main()
