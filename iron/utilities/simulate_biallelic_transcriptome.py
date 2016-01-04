#!/usr/bin/python
import argparse, sys, re, random, os, gzip
from multiprocessing import Pool, Lock, cpu_count, Queue
from subprocess import Popen, PIPE
import SimulationBasics
from VCFBasics import VCF
from SequenceBasics import read_fasta_into_hash
from TranscriptomeBasics import Transcriptome
from GenePredBasics import GenePredEntry
from RangeBasics import Bed,Locus,Loci
from FASTQPrecomputedProfileBasics import default_illumina_1_9 as default_illumina, default_pacbio_ccs95, default_pacbio_subreads
import FASTQBasics

write_lock = Lock()
gcounter = 0
shand1 = None #handles for writing out short reads
shand2 = None
emissions_reports = []

def main():
  parser = argparse.ArgumentParser(description="Create a simulated RNA-seq dataset")
  group0 = parser.add_mutually_exclusive_group(required=True)
  group0.add_argument('--load_biallelic_transcriptome',help="SERIALIZED BIALLELIC TRANSCRIOTOME EMITTER FILE to load up and use instead of all other file inputs")
  group0.add_argument('--inputs',nargs=3,help="<reference_genome> <phased_VCF> <transcripts_genepred>")
  #parser.add_argument('reference_genome',help="The reference genome.")
  #parser.add_argument('phased_VCF',help="A phased VCF file.  If you are simulating the genomes that step can make on of these for you.")
  #parser.add_argument('transcripts_genepred',help="A genepred file describing the transcripts.  Each transcript name must be unique.")
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--uniform_expression',action='store_true',help="Uniform distribution of transcript expression")
  group.add_argument('--isoform_expression',help="The transcript expression in TSV format <Transcript name> tab <Expression>")
  group.add_argument('--cufflinks_isoform_expression',help="The expression of the isoforms or - for a uniform distribution of transcript expression")
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--ASE_identical',type=float,help="The ASE for the transcriptome, every isoform will have the same allele preference.")
  group2.add_argument('--ASE_isoform_random',action='store_true',help="The ASE will be random for every isoform.")
  group2.add_argument('--ASE_locus_random',action='store_true',help="The ASE will be randomly assigned for each locus")
  parser.add_argument('--short_read_count',type=int,default=10000,help="INT number of short reads")
  parser.add_argument('--short_read_length',type=int,default=101,help="INT length of the short reads")
  parser.add_argument('--long_read_ccs_count',type=int,default=4000,help="INT default number of long reads")
  parser.add_argument('--long_read_subread_count',type=int,default=4000,help="INT default number of long reads")
  parser.add_argument('--no_errors',action='store_true',help="Do not simulate errors in reads")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Number of threads defaults to cpu_count()")
  parser.add_argument('--locus_by_gene_name',action='store_true',help="Faster than the complete calculation for overlapping loci.")
  parser.add_argument('--seed',type=int,help="seed to make transcriptome and rho creation deterministic.  Reads are still random, its just the transcriptome and rho that become determinisitic.")
  group3 = parser.add_mutually_exclusive_group(required=True)
  group3.add_argument('--output',help="Directory name for output")
  group3.add_argument('--save_biallelic_transcriptome',help="FILENAME output the biallelic transcriptome used to this file and then exit")
  parser.add_argument('--starting_read_multiplier',type=int,default=0,help="Used if outputting different reads from object, and you want them number differently give each different set values 0, 1, 2, etc...")
  args = parser.parse_args()
  fq_prof_illumina = None
  fq_prof_pacbio_ccs95 = None
  fq_prof_pacbio_subreads = None
  if not args.no_errors:
    fq_prof_illumina = default_illumina()
    fq_prof_pacbio_ccs95 = default_pacbio_ccs95()
    fq_prof_pacbio_subreads = default_pacbio_subreads()

  rbe = None
  if not args.load_biallelic_transcriptome:
    # we need to establish the emitter based on some known data
    rbe = load_from_inputs(args)
  
  else:
    rbe = SimulationBasics.RandomBiallelicTranscriptomeEmitter()
    inf = open(args.load_biallelic_transcriptome)
    sline = inf.readline().rstrip()
    inf.close()
    rbe.read_serialized(sline)
    
  if args.save_biallelic_transcriptome:
    ofser = open(args.save_biallelic_transcriptome,'w')
    ofser.write(rbe.get_serialized())
    ofser.close()
    return #exiting here
  # Lets prepare to output now
  args.output = args.output.rstrip('/')
  if not os.path.exists(args.output):
    os.makedirs(args.output)
  ofser = open(args.output+"/RandomBiallelicTranscriptomeEmitter.serialized",'w')
  ofser.write(rbe.get_serialized())
  ofser.close()
  rbe.set_gaussian_fragmentation_default_hiseq()
  #rbe_ser = rbe.get_serialized()
  sys.stderr.write("Sequencing short reads\n")
  global shand1
  shand1 = gzip.open(args.output+"/SR_1.fq.gz",'wb')
  global shand2
  shand2 = gzip.open(args.output+"/SR_2.fq.gz",'wb')
  z = 0
  buffer_full_size = 5000
  buffer = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for i in range(args.short_read_count*args.starting_read_multiplier,args.short_read_count*(args.starting_read_multiplier+1)):
    z = i+1
    buffer.append(z)
    if buffer_full_size <= len(buffer):
      vals = buffer[:]
      buffer = []
      if args.threads > 1:
        p.apply_async(process_short_read_buffer,args=(rbe,vals,args),callback=write_short_reads)
      else:
        oval = process_short_read_buffer(rbe,vals,args)
        write_short_reads(oval)
  if len(buffer) > 0:
    vals = buffer[:]
    buffer = []
    if args.threads > 1:
      p.apply_async(process_short_read_buffer,args=(rbe,vals,args),callback=write_short_reads)
    else:
      oval = process_short_read_buffer(rbe,vals,args)
      write_short_reads(oval)
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\nFinished sequencing short reads\n")
  shand1.close()
  shand2.close()
  global emissions_reports
  for i in range(0,len(emissions_reports)): emissions_reports[i]= emissions_reports[i].get()
  sr_report = combine_reports(emissions_reports)
  rbe.emissions_report = {} # initialize so we don't accidentally overwrite 
  # Now lets print out some of the emission details
  of = open(args.output+"/SR_report.txt",'w')
  for name in sorted(rbe.name2locus.keys()):
    express = 1
    if rbe.transcriptome1.expression:
      express = rbe.transcriptome1.expression.get_expression(name)
    if name in sr_report:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(sr_report[name][0])+"\t"+str(sr_report[name][1])+"\n")
    else:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(0)+"\t"+str(0)+"\n")
  of.close()


  rbe.emissions_report = {}
  emissions_reports = []
  # Now lets create the long read set
  rbe.set_gaussian_fragmentation_default_pacbio()
  #rbe_ser = rbe.get_serialized()
  sys.stderr.write("Sequencing long ccs reads\n")
  shand1 = gzip.open(args.output+"/LR_ccs95.fq.gz",'wb')
  buffer_full_size = 500
  buffer = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for i in range(args.starting_read_multiplier*args.long_read_ccs_count,(args.starting_read_multiplier+1)*args.long_read_ccs_count):
    z = i+1
    buffer.append(z)
    if buffer_full_size <= len(buffer):
      vals = buffer[:]
      buffer = []
      if args.threads > 1:
        p.apply_async(process_long_ccs_read_buffer,args=(rbe,vals,args),callback=write_long_reads)
      else:
        oval = process_long_ccs_read_buffer(rbe,vals,args)
        write_long_reads(oval)
  if len(buffer) > 0:
    vals = buffer[:]
    buffer = []
    if args.threads > 1:
      p.apply_async(process_long_ccs_read_buffer,args=(rbe,vals,args),callback=write_long_reads)
    else:
      oval = process_long_ccs_read_buffer(rbe,vals,args)
      write_long_reads(oval)
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\nFinished sequencing long reads\n")
  shand1.close()
  for i in range(0,len(emissions_reports)): emissions_reports[i]= emissions_reports[i].get()
  lr_ccs_report = combine_reports(emissions_reports)
  rbe.emissions_report = {} # initialize so we don't accidentally overwrite 
  # Now lets print out some of the emission details
  of = open(args.output+"/LR_ccs95_report.txt",'w')
  for name in sorted(rbe.name2locus.keys()):
    express = 1
    if rbe.transcriptome1.expression:
      express = rbe.transcriptome1.expression.get_expression(name)
    if name in lr_ccs_report:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(lr_ccs_report[name][0])+"\t"+str(lr_ccs_report[name][1])+"\n")
    else:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(0)+"\t"+str(0)+"\n")
  of.close()

  rbe.emissions_report = {}
  emissions_reports = []
  # Now lets create the long subread read set
  rbe.set_gaussian_fragmentation_default_pacbio()
  #rbe_ser = rbe.get_serialized()
  sys.stderr.write("Sequencing long subreads\n")
  shand1 = gzip.open(args.output+"/LR_subreads.fq.gz",'wb')
  buffer_full_size = 500
  buffer = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for i in range(args.long_read_subread_count*args.starting_read_multiplier,(args.starting_read_multiplier+1)*args.long_read_subread_count):
    z = i+1
    buffer.append(z)
    if buffer_full_size <= len(buffer):
      vals = buffer[:]
      buffer = []
      if args.threads > 1:
        p.apply_async(process_long_sub_read_buffer,args=(rbe,vals,args),callback=write_long_reads)
      else:
        oval = process_long_sub_read_buffer(rbe,vals,args)
        write_long_reads(oval)
  if len(buffer) > 0:
    vals = buffer[:]
    buffer = []
    if args.threads > 1:
      p.apply_async(process_long_sub_read_buffer,args=(rbe,vals,args),callback=write_long_reads)
    else:
      oval = process_long_sub_read_buffer(rbe,vals,args)
      write_long_reads(oval)
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\nFinished sequencing long reads\n")
  shand1.close()
  for i in range(0,len(emissions_reports)): emissions_reports[i]= emissions_reports[i].get()
  lr_sub_report = combine_reports(emissions_reports)
  rbe.emissions_report = {} # initialize so we don't accidentally overwrite 
  # Now lets print out some of the emission details
  of = open(args.output+"/LR_subreads_report.txt",'w')
  for name in sorted(rbe.name2locus.keys()):
    express = 1
    if rbe.transcriptome1.expression:
      express = rbe.transcriptome1.expression.get_expression(name)
    if name in lr_sub_report:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(lr_sub_report[name][0])+"\t"+str(lr_sub_report[name][1])+"\n")
    else:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(0)+"\t"+str(0)+"\n")
  of.close()

  combo_report = combine_reports([sr_report,lr_ccs_report,lr_sub_report])
  of = open(args.output+"/LR_SR_combo_report.txt",'w')
  for name in sorted(rbe.name2locus.keys()):
    express = 1
    if rbe.transcriptome1.expression:
      express = rbe.transcriptome1.expression.get_expression(name)
    if name in combo_report:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(combo_report[name][0])+"\t"+str(combo_report[name][1])+"\n")
    else:
      of.write(name +"\t"+rbe.gene_names[name]+"\t"+str(rbe.name2locus[name])+"\t"+str(express)+"\t"+str(rbe.transcriptome1_rho[name])+"\t"+str(0)+"\t"+str(0)+"\n")
  of.close()


def combine_reports(reports):
  c = {}
  for report in reports:
    for e in report:
      if e not in c: c[e] = [0,0]
      c[e][0] += report[e][0]
      c[e][1] += report[e][1]
  return c

def process_short_read_buffer(rbe,buffer,args):
  #rbe = SimulationBasics.RandomBiallelicTranscriptomeEmitter()
  #rbe.read_serialized(rbe_ser)
  fq_prof_illumina = default_illumina()
  read1 = ''
  read2 = ''
  zend = 0 
  for z in buffer:
    [name,l1,r1] = rbe.emit_paired_short_read(args.short_read_length)
    zend = z
    read1 += "@SRSIM"+str(z)+"\n"
    if args.no_errors:
      read1 += l1+"\n"
      read1 += "+\n"
      read1 += len(l1)*'I'+"\n"
    else:
      l1perm = fq_prof_illumina.create_fastq_and_permute_sequence(l1)
      read1 += l1perm['seq']+"\n"
      read1 += "+\n"
      read1 += l1perm['qual']+"\n"
    read2 += "@SRSIM"+str(z)+"\n"
    if args.no_errors:
      read2 += r1+"\n"
      read2 += "+\n"
      read2 += len(r1)*'I'+"\n"
    else:
      r1perm = fq_prof_illumina.create_fastq_and_permute_sequence(r1)
      read2 += r1perm['seq']+"\n"
      read2 += "+\n"
      read2 += r1perm['qual']+"\n"
  return [read1,read2,len(buffer),rbe.emissions_report]

def process_long_ccs_read_buffer(rbe,buffer,args):
  #rbe = SimulationBasics.RandomBiallelicTranscriptomeEmitter()
  #rbe.read_serialized(rbe_ser)
  fq_prof_pacbio_ccs95 = default_pacbio_ccs95()
  read1 = ''
  zend = 0 
  for z in buffer:
    [name,seq] = rbe.emit_long_read()
    g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
    zend = z
    read1 += "@"+g+"\n"
    if args.no_errors:
      read1 += seq+"\n"
      read1 += "+\n"
      read1 += len(seq)*'I'+"\n"
    else:
      seqperm = fq_prof_pacbio_ccs95.create_fastq_and_permute_sequence(seq)
      read1 += seqperm['seq']+"\n"
      read1 += "+\n"
      read1 += seqperm['qual']+"\n"
  return [read1,len(buffer),rbe.emissions_report]

def process_long_sub_read_buffer(rbe,buffer,args):
  #rbe = SimulationBasics.RandomBiallelicTranscriptomeEmitter()
  #rbe.read_serialized(rbe_ser)
  fq_prof_pacbio_subreads = default_pacbio_subreads()
  read1 = ''
  zend = 0 
  for z in buffer:
    [name,seq] = rbe.emit_long_read()
    g = 'm150102_010102_11112_c111111111111111112_s1_p0/'+str(z)+'/0_'+str(len(seq)-1)
    zend = z
    read1 += "@"+g+"\n"
    if args.no_errors:
      read1 += seq+"\n"
      read1 += "+\n"
      read1 += len(seq)*'I'+"\n"
    else:
      seqperm = fq_prof_pacbio_subreads.create_fastq_and_permute_sequence(seq)
      read1 += seqperm['seq']+"\n"
      read1 += "+\n"
      read1 += seqperm['qual']+"\n"
  return [read1,len(buffer),rbe.emissions_report]

def write_short_reads(vals):
  [read1,read2,zsize,emissions_report] = vals
  global write_lock
  global gcounter
  write_lock.acquire()
  global shand1
  global shand2
  global emissions_reports
  gcounter += zsize
  sys.stderr.write(str(gcounter)+"\r")
  shand1.write(read1)
  shand2.write(read2)
  eq = Queue()
  eq.put(emissions_report)
  emissions_reports.append(eq)
  write_lock.release()
  return

def write_long_reads(vals):
  [read1,zsize,emissions_report] = vals
  global write_lock
  global gcounter
  write_lock.acquire()
  global shand1
  global emissions_reports
  gcounter += zsize
  sys.stderr.write(str(gcounter)+"\r")
  shand1.write(read1)
  eq = Queue()
  eq.put(emissions_report)
  emissions_reports.append(eq)
  write_lock.release()
  return

# Pre:
# Take the allele info for one chromosome,
# Take one chromosome sequence string
# Take the left or right 0 or 1 position or the phased allele
# Take the chromosome name 
# Post:
# Reutrn a number of changes made, the chromosome name, and the chromosome sequence
def adjust_reference_genome(ainfo,refchrom,lrpos,chrom_name):
  reflist = list(refchrom)
  counter = 0
  for pos in sorted(ainfo):
    if reflist[pos-1].upper() != ainfo[pos][lrpos].upper():
      counter += 1
    reflist[pos-1] = ainfo[pos][lrpos]
  return [counter,chrom_name,''.join(reflist)]


def get_loci(transcripts_genepred):
  loci = Loci()
  loci.verbose= True
  with open(transcripts_genepred) as inf:
    for line in inf:
      if line[0]=='#': continue
      gpd = GenePredEntry(line.rstrip())
      rng = Bed(gpd.value('chrom'),gpd.value('txStart'),gpd.value('txEnd'))
      rng.set_payload(gpd.value('name'))
      loc1 = Locus()
      loc1.add_member(rng)
      loci.add_locus(loc1)
  sys.stderr.write("Organizing genepred data into overlapping loci\n")
  sys.stderr.write("Started with "+str(len(loci.loci))+" loci\n")
  loci.update_loci()
  sys.stderr.write("Ended with "+str(len(loci.loci))+" loci\n")

  m = 0
  locus2name = {}
  name2locus = {}
  for locus in loci.loci:
    m+=1
    for member in locus.members:
      name = member.get_payload()
      if m not in locus2name:  locus2name[m] = set()
      locus2name[m].add(name)
      name2locus[name] = m
  return [locus2name,name2locus]


def load_from_inputs(args):
  #Read in the VCF file
  sys.stderr.write("Reading in the VCF file\n")
  alleles = {}
  #with open(args.phased_VCF) as inf:
  with open(args.inputs[1]) as inf:
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

  sys.stderr.write("Reading in the reference genome\n")
  #ref = read_fasta_into_hash(args.reference_genome)
  ref = read_fasta_into_hash(args.inputs[0])
  res1 = []
  res2 = []
  p = None
  sys.stderr.write("Introducing VCF changes to reference sequences\n")
  # Pretty memory intesnive to so don't go with all possible threads
  if args.threads > 1: p = Pool(processes=max(1,int(args.threads/4)))
  for chrom in ref:
    # handle the case where there is no allele information
    if chrom not in alleles:
      r1q = Queue()
      r1q.put([0,chrom,ref[chrom]])
      res1.append(r1q)
      r2q = Queue()
      r2q.put([0,chrom,ref[chrom]])
      res2.append(r2q)
    elif args.threads > 1:
      res1.append(p.apply_async(adjust_reference_genome,args=(alleles[chrom],ref[chrom],0,chrom)))
      res2.append(p.apply_async(adjust_reference_genome,args=(alleles[chrom],ref[chrom],1,chrom)))
    else:
      r1q = Queue()
      r1q.put(adjust_reference_genome(alleles[chrom],ref[chrom],0,chrom))
      res1.append(r1q)
      r2q = Queue()
      r2q.put(adjust_reference_genome(alleles[chrom],ref[chrom],1,chrom))
      res2.append(r2q)
  if args.threads > 1:
    p.close()
    p.join()

  # now we can fill reference 1 with all our new sequences
  ref1 = {} 
  c1 = 0
  for i in range(0,len(res1)):
    res = res1[i].get()
    c1 += res[0]
    ref1[res[1]]=res[2]

  # now we can fill reference 2 with all our new sequences
  ref2 = {} 
  c2 = 0
  for i in range(0,len(res2)):
    res = res2[i].get()
    c2 += res[0]
    ref2[res[1]]=res[2]
  sys.stderr.write("Made "+str(c1)+"|"+str(c2)+" changes to the reference\n")

  # Now ref1 and ref2 have are the diploid sources of the transcriptome
  gpdnames = {}
  txn1 = Transcriptome()
  txn2 = Transcriptome()
  txn1.set_reference_genome_dictionary(ref1)
  txn2.set_reference_genome_dictionary(ref2)
  #with open(args.transcripts_genepred) as inf:
  with open(args.inputs[2]) as inf:
    for line in inf:
      if line[0]=='#': continue
      txn1.add_genepred_line(line.rstrip())
      txn2.add_genepred_line(line.rstrip())
      gpd = GenePredEntry(line.rstrip())
      gpdnames[gpd.value('name')] = gpd.value('gene_name')
  # The transcriptomes are set but we dont' really need the references anymore
  # Empty our big memory things
  txn1.ref_hash = None
  txn2.ref_hash = None
  for chrom in ref1.keys():  del ref1[chrom]
  for chrom in ref2.keys():  del ref2[chrom]
  for chrom in ref.keys():  del ref[chrom]

  if not args.locus_by_gene_name:
    #[locus2name,name2locus] = get_loci(args.transcripts_genepred)
    [locus2name,name2locus] = get_loci(args.inputs[2])
  else: # set locus by gene name
    sys.stderr.write("Organizing loci by gene name\n")
    locus2name = {}
    name2locus = {}
    numname = {}
    m = 0
    for name in sorted(gpdnames): 
      gene = gpdnames[name]
      if gene not in numname:
        m+=1
        numname[gene] = m
      num = numname[gene]
      if num not in locus2name:
        locus2name[num] = set()
      locus2name[num].add(name)
      name2locus[name] = num
    sys.stderr.write("Ended with "+str(len(locus2name.keys()))+" loci\n")

  if args.isoform_expression:
    sys.stderr.write("Reading expression from a TSV\n")
    with open(args.isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        txn1.add_expression(f[0],float(f[1]))
        txn2.add_expression(f[0],float(f[1]))
  elif args.cufflinks_isoform_expression:
    sys.stderr.write("Using cufflinks expression\n")
    cuffz = 0
    with open(args.cufflinks_isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        cuffz +=1
        sys.stderr.write(str(cuffz)+" cufflinks entries processed\r")
        f = line.rstrip().split("\t")
        txn1.add_expression_no_update(f[0],float(f[9]))
        txn2.add_expression_no_update(f[0],float(f[9]))
    txn1.update_expression()
    txn2.update_expression()
    sys.stderr.write("\n")
  elif args.uniform_expression:
    sys.stderr.write("Using uniform expression model\n")
  else:
    sys.stderr.write("Warning isoform expression not sepcified, using uniform expression model.\n")
  # Now we have the transcriptomes set
  rhos = {} # The ASE of allele 1 (the left side)
  randos = {}
  if args.seed:
    random.seed(args.seed)
  for z in locus2name: randos[z] = random.random()
  sys.stderr.write("Setting rho for each transcript\n")
  # Lets set rho for ASE for each transcript
  for tname in sorted(txn1.transcripts):
    if args.ASE_identical or args.ASE_identical == 0:
      rhos[tname] = float(args.ASE_identical)
    elif args.ASE_isoform_random:
      rhos[tname] = random.random()
    else: # we must be on locus random
      rhos[tname] = randos[name2locus[tname]]
  #Now our dataset is set up
  rbe = SimulationBasics.RandomBiallelicTranscriptomeEmitter(txn1,txn2)
  rbe.gene_names = gpdnames
  rbe.name2locus = name2locus
  rbe.set_transcriptome1_rho(rhos)
  return rbe

if __name__=="__main__":
  main()
