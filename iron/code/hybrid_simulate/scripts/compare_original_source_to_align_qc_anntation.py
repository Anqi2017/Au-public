#!/usr/bin/python
import sys, argparse, gzip

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('original_source',help="Generated from --output_original_source in hybrid_simulate emit")
  parser.add_argument('align_qc_annotation',help="Stored in annotbest.txt.gz")
  parser.add_argument('--original_source_fields',default=[0,1,2],nargs=3,type=int,help="zero index fields <read> <gene> <tx>")
  parser.add_argument('--align_qc_annotation_fields',default=[1,2,3,4],nargs=4,type=int,help="zero index fields <read> <gene> <tx> <type>")
  args = parser.parse_args()
  
  inf1 = None
  if args.original_source[-3:] == '.gz':
    inf1 = gzip.open(args.original_source)
  else:
    inf1 = open(args.original_source)
  reads = {}
  for line in inf1:
    f = line.rstrip().split("\t")
    read = f[args.original_source_fields[0]]
    gene = f[args.original_source_fields[1]]
    tx = f[args.original_source_fields[2]]
    reads[f[0]] = [gene,tx]
  inf1.close()
  inf2 = None
  if args.align_qc_annotation[-3:] == '.gz':
    inf2 = gzip.open(args.align_qc_annotation)
  else:
    inf2 = open(args.align_qc_annotation)
  annot = {}
  for line in inf2:
    f = line.rstrip().split("\t")
    read = f[args.align_qc_annotation_fields[0]]
    gene = f[args.align_qc_annotation_fields[1]]
    tx = f[args.align_qc_annotation_fields[2]]
    type = f[args.align_qc_annotation_fields[3]]
    annot[read] = [gene, tx, type]
  inf2.close()
  unannotated = []
  txcorrect_full = 0
  txcorrect_partial = 0
  txwrong = 0
  genecorrect = 0
  genewrong = 0
  annotated = 0
  total = 0
  for read in reads:
    total += 1
    if read not in annot:
      unannotated.append(read)
      continue
    annotated += 1
    (agene,atx,atype) = annot[read]
    if atx == reads[read][1] and (atype == 'full' or atype=='Full'): txcorrect_full+=1
    elif atx == reads[read][1] and (atype == 'partial' or atype=='Partial'): txcorrect_partial+=1
    else: txwrong += 1
    if agene == reads[read][0]: genecorrect += 1
    else: genewrong += 1
  print "Annotated\t"+str(annotated)+"/"+str(total)+"\t"+str(float(annotated)/float(total))
  print "CorrectTx\t"+str(txcorrect_full+txcorrect_partial)+"/"+str(total)+"\t"+str(float(txcorrect_full+txcorrect_partial)/float(total))
  print "FullTx\t"+str(txcorrect_full)
  print "PartialTx\t"+str(txcorrect_partial)
  print "CorrectGene\t"+str(genecorrect)+"/"+str(total)+"\t"+str(float(genecorrect)/float(total))
if __name__=="__main__":
  main()
