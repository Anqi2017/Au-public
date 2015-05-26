#!/usr/bin/python

import sys, os, inspect, argparse

#pre: psl filename, an optional integer for smoothing size
#     The smoothing is done to combine adjacent "exons" within that size definition to get rid of small deletions or indels
#post: a genepred file with target sequence information.  not query

#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

import PSLBasics, GenePredBasics


def main():
  parser = argparse.ArgumentParser(description="Convert a psl file into a target formated genepred file.")
  parser.add_argument('--fill_gaps',type=int,default=0,help="Close gaps this size or smaller.")
  parser.add_argument('input_name',help="Input PSL file, use - to indicate STDIN.")
  args = parser.parse_args()
  
  pslfilehandle = sys.stdin
  if args.input_name != '-':
    pslfilehandle = open(args.input_name)
  with pslfilehandle as infile:
    for line in infile:
      psl_entry = PSLBasics.line_to_entry(line)
      genepred_line = PSLBasics.convert_entry_to_genepred_line(psl_entry)
      if args.fill_gaps > 0:
        genepred_entry = GenePredBasics.line_to_entry(genepred_line)
        genepred_entry2 = GenePredBasics.smooth_gaps(genepred_entry,args.fill_gaps)
        genepred_line = GenePredBasics.entry_to_line(genepred_entry2)
      print genepred_line

if __name__=="__main__":
  main()
