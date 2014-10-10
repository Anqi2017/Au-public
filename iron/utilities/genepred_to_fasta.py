#!/usr/bin/python
import sys, os, inspect

#pre: a genepred file, a reference fasta file, an output fasta, optionally 'directionless' will reverse compliment entries on the negative strand
#post: writes to the output fasta the sequence of the transcripts in the genepred
#      produces transcripts in the proper orientation according to their direction in the genepred
#modifies fileIO

#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)


import genepred_basics


def main():

  if len(sys.argv) < 4:
    print sys.argv[0] + " <genepred> <genome fasta> <output fasta> <(optional) 'directionless' for no RC>"
    sys.exit()
  genepred_filename = sys.argv[1]
  genome_filename = sys.argv[2]
  output_filename = sys.argv[3]
  dodirectionless = 0
  if len(sys.argv) == 5:
    dodirectionless = 1
  if dodirectionless == 1:
    genepred_basics.write_genepred_to_fasta_directionless(genepred_filename,genome_filename,output_filename)
  else:
    genepred_basics.write_genepred_to_fasta(genepred_filename,genome_filename,output_filename)

main()
