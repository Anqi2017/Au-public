#!/usr/bin/python

import sys, os, inspect, re

#pre: fasta for 
#     1.  CCS reads called with 95% accuracy 
#     2.  CCS reads called with 90% accuracy
#     3.  Subreads called with 75% accuracy
#post: a fasta excluding reads in the ccs 95 set
#      with ccs reads with accuracy from 90-95 and the 
#      longest subreads with 75 or better

#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)
import sequence_basics

def main():
  if len(sys.argv) < 4:
    print sys.argv[0] + ' <ccs 95 fasta> <ccs 90 fasta> <subread 75 fasta>'
    sys.exit()
  ccs95 = sequence_basics.read_fasta_into_array(sys.argv[1])

  #get the names and base names from ccs95
  ccs95basenames = set()
  prog = re.compile('(^[^\/]+\/\d+)\/')
  for entry in ccs95:
    m = prog.match(entry['name'])
    if not m: 
      print 'trouble parsing name'
      return
    basename = m.group(1)
    ccs95basenames.add(basename)
  ccs95 = []

  ccs90 = sequence_basics.read_fasta_into_array(sys.argv[2])
  ccs90to95basenames = set()
  for entry in ccs90:
    m = prog.match(entry['name'])
    if not m:
      print 'trouble parsing name'
      return
    basename = m.group(1)
    if basename in ccs95basenames: continue
    print '>'+entry['name']
    print entry['seq']
    ccs90to95basenames.add(basename)
  ccs90 = []

  sub75 = sequence_basics.read_fasta_into_array(sys.argv[3])
  sub75lengths = {}
  for entry in sub75:
    m = prog.match(entry['name'])
    if not m:
      print 'trouble parsing name'
      return
    basename = m.group(1)
    if basename in ccs95basenames: continue
    if basename in ccs90to95basenames: continue
    if basename not in sub75lengths: 
      sub75lengths[basename] = len(entry['seq'])
  printsub75 = {}
  for entry in sub75:
    m = prog.match(entry['name'])
    if not m:
      print 'trouble parsing name'
      return
    basename = m.group(1)
    if basename in ccs95basenames: continue
    if basename in ccs90to95basenames: continue
    if len(entry['seq']) == sub75lengths[basename]:
      printsub75[basename] = entry
  for basename in printsub75:
    entry = printsub75[basename]
    print '>'+entry['name']
    print entry['seq']
  sys.stderr.write('"Dear Benjamin: Everything thing finished nicely and its all going to be okay now." - STDERR')
    
main()
