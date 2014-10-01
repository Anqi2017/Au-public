#!/usr/bin/python
import sys
import re

################
#  Jason Weirather 20140317
#  Get chromosome names as the first non-whitespace characters
#  Pre:  fasta file
#  Post:  write list of chromosomes to standard output
#  Modifies:  standard output
################

if(len(sys.argv) < 2):
  print 'split_genome_fasta.py <INPUT FILE>'
  sys.exit()

with open(sys.argv[1]) as fp:
  for line in fp:
    if line.startswith('>'):
      p = re.compile('>\s*(\S+)')
      name = p.match(line).group(1)
      print name
