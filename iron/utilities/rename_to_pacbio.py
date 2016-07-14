#!/usr/bin/python
import sys,argparse
from SequenceBasics import FastaHandleReader, FastqHandleReader

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('input',help="Use - for STDIN")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--fasta',action='store_true')
  group.add_argument('--fastq',action='store_true')
  group.add_argument('--gpd',action='store_true')
  parser.add_argument('--output_table',help='save coversion to file')
  parser.add_argument('-o','--output')
  args = parser.parse_args()

  if args.input=='-': args.input = sys.stdin
  else: args.input= open(args.input)

  if args.output: args.output = open(args.output,'w')
  if args.output_table:  args.output_table= open(args.output_table,'w')
  else: args.output = sys.stdout
  if args.gpd:
    z = 0
    for line in args.input:
        f = line.rstrip().split("\t")
        z+=1
        name = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
        if args.output_table: args.output_table.write(f[0]+"\t"+name+"\n")
        f[0] = name
        f[1] = name
        args.output.write("\t".join(f)+"\n")
    args.output.close()
    if args.output_table:
      args.output_table.close()
    return


  if args.fasta:
    args.input = FastaHandleReader(args.input)
  elif args.fastq:
    args.input = FastqHandleReader(args.input)
  z = 0
  while True:
    e = args.input.read_entry()
    if not e: break
    z+=1
    name = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
    if args.fastq:
      args.output.write( '@'+name+"\n"+ e['seq']+"\n"+ '+'+e['qual']+"\n")
    elif args.fasta:
      args.output.write('>'+name+"\n"+e['seq']+"\n")
    if args.output_table: args.output_table.write(e['name']+"\t"+name+"\n")
  args.output.close()
  if args.output_table: args.output_table.close()
if __name__=="__main__":
  main()
