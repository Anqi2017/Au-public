#!/usr/bin/python
import sys, argparse, re, os

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('inputs',nargs='+')
  parser.add_argument('--prefix',help="add to begining")
  args = parser.parse_args()
  
  for input in args.inputs:
    m = re.match('^(.*)\.sorted\.bam$',input)
    if not m: 
      sys.stderr.write("ERROR unexpected input file type "+input+"\n")
      sys.exit()
    if os.path.exists(m.group(1)+'.alignqc-gencode.xhtml'): 
      sys.stderr.write("WARNING already finished "+m.group(1)+" skipping\n")
      continue
    cmd = 'alignqc analyze --context_error_stopping_point 10000 --alignment_error_max_length 1000000 '
    cmd += ' -r /Shared/Au/jason/Reference/UCSC/Human/hg20_GRCh38_dec2013/Genome/genome.fa '
    cmd += ' -a /Shared/Au/jason/Reference/UCSC/Human/hg20_GRCh38_dec2013/Genes/genePred/gencode.v24.annotation.20160712.gpd '
    cmd += input
    cmd += ' -o '+m.group(1)+'.alignqc-gencode.xhtml '
    cmd += ' --portable_output '+m.group(1)+'.alignqc-gencode.portable.xhtml '
    cmd += ' --output_folder '+m.group(1)+'.alignqc-gencode '
    cmd += ' --threads 4 '
    cmd += ' --tempdir /localscratch/Users/weirathe '
    if args.prefix:
      cmd = args.prefix + ' ' +cmd
    print cmd

if __name__=="__main__":
  main()
