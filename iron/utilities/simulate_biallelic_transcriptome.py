#!/usr/bin/python
import argparse, sys, re, random, os, gzip
from subprocess import Popen, PIPE
import SimulationBasics
from VCFBasics import VCF
from SequenceBasics import read_fasta_into_hash
from TranscriptomeBasics import Transcriptome

def main():
  parser = argparse.ArgumentParser(description="Create a simulated RNA-seq dataset")
  parser.add_argument('reference_genome',help="The reference genome.")
  parser.add_argument('phased_VCF',help="A phased VCF file.  If you are simulating the genomes that step can make on of these for you.")
  parser.add_argument('transcripts_genepred',help="A genepred file describing the transcripts.  Each transcript name must be unique.")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--uniform_expression',action='store_true',help="Uniform distribution of transcript expression")
  group.add_argument('--isoform_expression',help="The transcript expression in TSV format <Transcript name> tab <Expression>")
  group.add_argument('--cufflinks_isoform_expression',help="The expression of the isoforms or - for a uniform distribution of transcript expression")
  parser.add_argument('ASE_expression',help="The ASE for the transcriptome, or float for a specifc ASE for all genes, or - for uniform distribution")
  parser.add_argument('output',help="Directory name for output")
  parser.add_argument('--short_read_count',type=int,default=10000,help="INT number of short reads")
  parser.add_argument('--short_read_length',type=int,default=101,help="INT length of the short reads")
  parser.add_argument('--long_read_count',type=int,default=4000,help="INT default number of long reads")
  args = parser.parse_args()
  args.output = args.output.rstrip('/')
  #Read in the VCF file
  alleles = {}
  with open(args.phased_VCF) as inf:
    for line in inf:
      vcf = VCF(line)
      if not vcf.is_snp(): continue
      g = vcf.get_phased_genotype()
      if not g: continue
      if vcf.value('chrom') not in alleles:
        alleles[vcf.value('chrom')] = {}
      if vcf.value('pos') in alleles[vcf.value('chrom')]:
        sys.stderr.write("WARNING: seeing the same position twice.\n"+line.rstrip()+"\n")
      alleles[vcf.value('chrom')][vcf.value('pos')] = g # set our left and right
  ref1 = read_fasta_into_hash(args.reference_genome)
  ref2 = ref1.copy() # copy the dictionary
  c1 = 0
  c2 = 0
  for chrom in alleles:
    rseq1 = list(ref1[chrom])
    rseq2 = list(ref2[chrom])
    ainfo = alleles[chrom]
    for pos in sorted(ainfo):
      if rseq1[pos-1].upper() != ainfo[pos][0].upper():
        c1 += 1
      rseq1[pos-1] = ainfo[pos][0]
      if rseq2[pos-1].upper() != ainfo[pos][1].upper():
        c2 += 1
      rseq2[pos-1] = ainfo[pos][1]
    ref1[chrom] = ''.join(rseq1)
    ref2[chrom] = ''.join(rseq2)
  sys.stderr.write("Made "+str(c1)+"|"+str(c2)+" changes to the reference\n")
  # Now ref1 and ref2 have are the diploid sources of the transcriptome
  txn1 = Transcriptome()
  txn2 = Transcriptome()
  txn1.set_reference_genome_dictionary(ref1)
  txn2.set_reference_genome_dictionary(ref2)
  with open(args.transcripts_genepred) as inf:
    for line in inf:
      if line[0]=='#': continue
      txn1.add_genepred_line(line.rstrip())
      txn2.add_genepred_line(line.rstrip())
  if args.isoform_expression:
    sys.stderr.write("Reading expression from a TSV\n")
    with open(args.isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        txn1.add_expression(f[0],float(f[1]))
        txn2.add_expression(f[0],float(f[1]))
  elif args.uniform_expression:
    sys.stderr.write("Using uniform expression model\n")
  elif args.cufflinks_isoform_expression:
    sys.stderr.write("Using cufflinks expression\n")
    with open(args.cufflinks_isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        txn1.add_expression(f[0],float(f[9]))
        txn2.add_expression(f[0],float(f[9]))
  # Now we have the transcriptomes set
  rhos = {} # The ASE of allele 1 (the left side)
  # Lets set rho for ASE for each transcript
  for tname in txn1.transcripts:
    m = re.match('([01]*\.*\d+)',args.ASE_expression)
    if args.ASE_expression == '-':
      rhos[tname] = random.random()
    elif m:
      rhos[tname] = float(m.group(1))
  #Now our dataset is set up
  rbe = SimulationBasics.RandomBiallelicTranscriptomeEmitter(txn1,txn2)
  rbe.set_transcriptome1_rho(rhos)
  rbe.set_gaussian_fragmentation_default_hiseq()
  # Lets prepare to output now
  if not os.path.exists(args.output):
    os.makedirs(args.output)
  sys.stderr.write("Sequencing short reads\n")
  of1 = gzip.open(args.output+"/SR_1.fq.gz",'wb')
  of2 = gzip.open(args.output+"/SR_2.fq.gz",'wb')
  for i in range(0,args.short_read_count):
    z = i+1
    if z %100==0: sys.stderr.write(str(z)+"\r")
    [name,l1,r1] = rbe.emit_paired_short_read(args.short_read_length)
    of1.write("@SRSIM"+str(z)+"\n")
    of1.write(l1+"\n")
    of1.write("+\n")
    of1.write(len(l1)*'I'+"\n")
    of2.write("@SRSIM"+str(z)+"\n")
    of2.write(r1+"\n")
    of2.write("+\n")
    of2.write(len(r1)*'I'+"\n")
  sys.stderr.write("\nFinished sequencing short reads\n")
  of1.close()
  of2.close()

  # Now lets print out some of the emission details
  of = open(args.output+"/SR_report.txt",'w')
  for name in sorted(rbe.emissions_report.keys()):
    of.write(name +"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(rbe.emissions_report[name][0])+"\t"+str(rbe.emissions_report[name][1])+"\n")
  of.close()
  rbe.emissions_report = {}

  # Now lets create the long read set
  rbe.set_gaussian_fragmentation_default_pacbio()
  sys.stderr.write("Sequencing long reads\n")
  of = gzip.open(args.output+"/LR.fq.gz",'wb')
  for i in range(0,args.long_read_count):
    z = i+1
    if z %100==0: sys.stderr.write(str(z)+"\r")
    [name,seq] = rbe.emit_long_read()
    g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
    of.write("@"+g+"\n")
    of.write(seq+"\n")
    of.write("+\n")
    of.write(len(seq)*'I'+"\n")
  of.close()
  sys.stderr.write("\nFinished sequencing long reads\n")
  # Now lets print out some of the emission details
  of = open(args.output+"/LR_report.txt",'w')
  for name in sorted(rbe.emissions_report.keys()):
    of.write(name +"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(rbe.emissions_report[name][0])+"\t"+str(rbe.emissions_report[name][1])+"\n")
  of.close()
  combo = {}
  with open(args.output+"/SR_report.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      [name,rho,left,right] = f
      if name not in combo:
        combo[name] = {}
        combo[name]['rho'] = rho
        combo[name]['left'] = 0
        combo[name]['right'] = 0
      combo[name]['left'] += int(left)
      combo[name]['right'] += int(right)
  with open(args.output+"/LR_report.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      [name,rho,left,right] = f
      if name not in combo:
        combo[name] = {}
        combo[name]['rho'] = rho
        combo[name]['left'] = 0
        combo[name]['right'] = 0
      combo[name]['left'] += int(left)
      combo[name]['right'] += int(right)
  of = open(args.output+"/LR_SR_combo_report.txt",'w')
  for name in sorted(combo):
    of.write(name+"\t"+combo[name]['rho']+"\t"+str(combo[name]['left'])+"\t"+str(combo[name]['right'])+"\n")
  of.close()
if __name__=="__main__":
  main()
