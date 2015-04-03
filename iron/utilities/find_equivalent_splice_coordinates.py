#!/usr/bin/python
import argparse, sys
from SequenceBasics import read_fasta_into_hash, rc

# Pre: A TSV file of splice site coordinates. Each 
#      line contains a pair of splice coordinates.
#      see the parser help for more information on 
#      format
#      A genome reference fasta file
# Post: A like-wise formatted TSV file that includes
#       a line for each equivalent splice.
#       If long_form switch is on, it will have an
#       extra column indicating whether the line is
#       an 'original' entry, or an 'alternate' equivalent
#       splice
#       Alternate only is an option also.
# Modifies: STDOUT

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('--only_output_alternates',action='store_true',help='When selected, the original coordiantes are not output, and only the alternates are output')
  parser.add_argument('--long_form', action='store_true',help="add an additional column to the beginning of the output indicating whether it is an original or alternate splice coordinate")
  parser.add_argument('GenomeFastaFile',nargs=1,help="FILENAME Fasta format file of the reference genome")
  parser.add_argument('SpliceSiteFile',nargs=1,help="FILENAME Splice Site file is in tsv format with <Left chrom> <Left coord (base-1)> <Left dir [+-]> <Right chrom> <Right coord (base-1)> <Right dir [+-]>\nWhere the coordinates indicate the base that is inside the exon proximal to the splice.  Direction indicates the transcription direction on the chromosome for that side of the splice.  For coordiantes 1-base means that the number 1 would be the first base of the sequence (makes sense to do it that way, right? :P)")
  of = sys.stdout
  args = parser.parse_args()
  golds = []
  with open(args.SpliceSiteFile[0]) as inf:
    for line in inf:
      f = line.rstrip().split()
      t = {}
      t['l'] = {}
      t['r'] = {}
      t['l']['chr'] = f[0]
      t['l']['coord'] = int(f[1])
      t['l']['dir'] = f[2]
      t['r']['chr'] = f[3]
      t['r']['coord'] = int(f[4])
      t['r']['dir'] = f[5]
      golds.append(t)

  ref = read_fasta_into_hash(args.GenomeFastaFile[0])
  lens = {}
  for chr in ref:
    lens[chr] = len(ref[chr])
  for g in golds:
    l_chrom = g['l']['chr']
    r_chrom = g['r']['chr']
    l_start = g['l']['coord']
    r_start = g['r']['coord']
    l_dir = g['l']['dir']
    r_dir = g['r']['dir']
    # print the main case
    if not args.only_output_alternates:
      startstring = ''
      if args.long_form: startstring = "original\t"
      of.write(startstring+l_chrom + "\t" + str(l_start) + "\t" + l_dir + "\t" + r_chrom + "\t" + str(r_start) + "\t" + r_dir+"\n")
    #check upstream left
    equivalent = 1
    l_base = l_start
    r_base = r_start
    while(equivalent == 1):
      left_bases = ''
      right_bases = ''
      if l_dir == '+':
        l_base -= 1
        if l_base < 1: break
        left_bases = str(ref[l_chrom][l_base])
      else:
        l_base += 1
        if l_base > lens[l_chrom]: break
        left_bases = rc(str(ref[l_chrom][l_base-2]))
      if r_dir == '+':
        r_base -= 1
        if r_base < 1: break
        right_bases = str(ref[r_chrom][r_base-1])
      else:
        r_base += 1
        if r_base > lens[r_chrom]: break
        right_bases = rc(str(ref[r_chrom][r_base-1]))
      if left_bases != right_bases: break
      startstring = ''
      if args.long_form: startstring = "alternate\t"
      of.write(startstring+l_chrom + "\t" + str(l_base) + "\t" + l_dir + "\t" + r_chrom + "\t" + str(r_base) + "\t" + r_dir+"\n")
    #check downstream left
    equivalent = 1
    l_base = l_start
    r_base = r_start
    while(equivalent == 1):
      left_bases = ''
      right_bases = ''
      if l_dir == '+':
        l_base += 1
        if l_base > lens[l_chrom]: break
        left_bases = str(ref[l_chrom][l_base-1])
      else:
        l_base -= 1
        if l_base < 1: break
        left_bases = rc(str(ref[l_chrom][l_base-1]))
      if r_dir == '+':
        r_base += 1
        if r_base > lens[r_chrom]: break
        right_bases = str(ref[r_chrom][r_base-2])
      else:
        r_base -= 1
        if r_base > lens[r_chrom]: break
        right_bases = rc(str(ref[r_chrom][r_base]))
      if left_bases != right_bases: break
      startstring = ''
      if args.long_form: startstring = "alternate\t"
      of.write(startstring+l_chrom + "\t" + str(l_base) + "\t" + l_dir + "\t" + r_chrom + "\t" + str(r_base) + "\t" + r_dir+"\n")

main()
