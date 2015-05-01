#!/usr/bin/python
import argparse, re, sys
from SequenceBasics import GenericFastaFileReader
# Get ccs reads first and then get the longest subreads that we don't have a subread for

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('CCS_file',nargs=1,help='FILENAME CCS file')
  parser.add_argument('subread_file',nargs=1,help='FILENAME subread file')
  args = parser.parse_args()
  ccs_gffr = GenericFastaFileReader(args.CCS_file[0])
  names = set()
  while True:
    entry = ccs_gffr.read_entry()
    if not entry: break
    m = re.match('^([^\/]+\/\d+)\/ccs', entry['name'])
    if not m: 
      sys.stderr.write("ERROR:  strange format for what should be a ccs read header\n")
      return
    names.add(m.group(1))
    print entry['original'].rstrip()
  ccs_gffr.close()
  sub_gffr = GenericFastaFileReader(args.subread_file[0])
  longest = {}
  longest_name = {}
  while True:
    entry = sub_gffr.read_entry()
    if not entry: break
    m = re.match('^([^\/]+\/\d+)\/', entry['name'])
    obs = m.group(1)
    if obs in names: continue #already have it from ccs
    if obs not in longest:
      longest[obs] = 0
    if len(entry['seq']) > longest[obs]:
      longest_name[obs] = entry['name']
  sub_gffr.close()
  sub_gffr = GenericFastaFileReader(args.subread_file[0])
  while True:
    entry = sub_gffr.read_entry()
    if not entry: break
    m = re.match('^([^\/]+\/\d+)\/', entry['name'])
    if m.group(1) in names: continue # can still skip the ones we have from ccs
    if entry['name'] == longest_name[m.group(1)]:
      print entry['original'].rstrip()
  sub_gffr.close()    

main()
