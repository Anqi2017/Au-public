#!/usr/bin/python

import sys, re, os, inspect, argparse

#pre: genepred filename, an optional integer for smoothing size
#     The smoothing is done to combine adjacent "exons" within that size definition to get rid of small deletions or indels
#post: a bed file where each line has chromosome, index 1 start, and index 1 end coordinates and the gene_name from each psl line is listed

#bring in the folder to the path for our modules
#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)
import GenePredBasics


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--closegap',nargs='?',help='INT close gaps less than or equal to this, if set.')
  parser.add_argument('--namefield',nargs='?',help='INT[1 or 2] use the first or second field. Default (1)')
  parser.add_argument('--noheader',help='do not print any tack header')
  parser.add_argument('--headername',nargs='?',help='STRING name for the track. Default is the gpd file name')
  parser.add_argument('--headerdescription',nargs='?',help='STRING description for the track.')
  parser.add_argument('--color',choices=['blue','green','orange','purple','red'])
  parser.add_argument('in_genepredfile',nargs=1,help='FILENAME input genepred, use - for STDIN')
  args = parser.parse_args()

  color = '0,0,0'

  if args.color:
    if args.color == 'blue':
      color = '67,162,202'
    elif args.color == 'green':
      color = '49,163,84'
    elif args.color == 'orange':
      color = '254,178,76'
    elif args.color == 'purple':
      color = '136,86,167'
    elif args.color == 'red':
      color = '240,59,32'

  genepredfilename = args.in_genepredfile[0]
  smooth_size = 0
  if args.closegap:
    smooth_size = int(args.closegap)
  namefield = 1
  if args.namefield:
    namefield = int(args.namefield)

  # set up the header if one is desired
  header = ''
  if not args.noheader:
    newname = 'longreads'
    m = re.search('([^\/]+)$',genepredfilename)
    if m:
      newname = m.group(1)
    newname = re.sub('[\s]+','_',newname)
    if args.headername:
      newname = args.headername
    elif genepredfilename == '-':
      newname = 'STDIN'
    header += "track\tname="+newname+"\t"
    description = newname+' GenePred Entries'
    if args.headerdescription:
       description = args.headerdescription
    header += 'description="'+description + '"'+"\t"
    header += 'itemRgb="On"'
    print header
  
  gpd_handle = sys.stdin
  if genepredfilename != '-': gpd_handle = open(genepredfilename)

  with gpd_handle as infile:
    for line in infile:
      if re.match('^#',line):
        continue
      genepred_entry = GenePredBasics.line_to_entry(line)
      if smooth_size > 0:
        genepred_entry = GenePredBasics.smooth_gaps(genepred_entry,smooth_size)
      exoncount = len(genepred_entry['exonEnds'])
      ostring  = genepred_entry['chrom'] + "\t" 
      ostring += str(genepred_entry['exonStarts'][0]) + "\t"
      ostring += str(genepred_entry['exonEnds'][exoncount-1]) + "\t"
      if namefield == 1:
        ostring += genepred_entry['gene_name'] + "\t"
      else: 
        ostring += genepred_entry['name']
      ostring += '1000' + "\t"
      ostring += genepred_entry['strand'] + "\t" 
      ostring += str(genepred_entry['exonStarts'][0]) + "\t"
      ostring += str(genepred_entry['exonEnds'][exoncount-1]) + "\t"      
      ostring += color+"\t"
      ostring += str(exoncount) + "\t"
      for i in range(0,exoncount):
        ostring += str(genepred_entry['exonEnds'][i]-genepred_entry['exonStarts'][i]) + ','
      ostring += "\t"
      for i in range(0,exoncount):
        ostring += str(genepred_entry['exonStarts'][i]-genepred_entry['exonStarts'][0])+','
      print ostring
      #for i in range(0,len(genepred_entry['exonStarts'])):
main()
