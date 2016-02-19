#!/usr/bin/python
import argparse, sys, re, random, os, gzip
from subprocess import Popen, PIPE
import SimulationBasics
from SequenceBasics import read_fasta_into_hash
from TranscriptomeBasics import Transcriptome
from FASTQPrecomputedProfileBasics import default_illumina_1_9 as default_illumina, default_pacbio_ccs95, default_pacbio_subreads
import FASTQBasics
from multiprocessing import Pool, Lock, cpu_count

left_handle = None
right_handle = None
long_handle = None
glock = Lock()
greport = {}

def main():
  parser = argparse.ArgumentParser(description="Create a simulated RNA-seq dataset")
  parser.add_argument('reference_genome',help="The reference genome.")
  parser.add_argument('transcripts_genepred',help="A genepred file describing the transcripts.  Each transcript name must be unique.")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--uniform_expression',action='store_true',help="Uniform distribution of transcript expression")
  group.add_argument('--isoform_expression',help="The transcript expression in TSV format <Transcript name> tab <Expression>")
  group.add_argument('--cufflinks_isoform_expression',help="The expression of the isoforms or - for a uniform distribution of transcript expression")
  group2 = parser.add_mutually_exclusive_group()
  group2.add_argument('--long_reads_only',action='store_true')
  group2.add_argument('--short_reads_only',action='store_true')
  group2.add_argument('--output',help="Directory name for output")
  parser.add_argument('--short_read_count',type=int,default=10000,help="INT number of short reads")
  parser.add_argument('--short_read_length',type=int,default=101,help="INT length of the short reads")
  parser.add_argument('--long_read_count',type=int,default=4000,help="INT default number of long reads")
  parser.add_argument('--no_errors',action='store_true')
  parser.add_argument('--threads',type=int,default=1)
  args = parser.parse_args()
  if args.output:
    args.output = args.output.rstrip('/')

  fq_prof_pacbio_ccs95 = None
  fq_prof_pacbio_subreads = None
  fq_prof_illumina = None
  if not args.no_errors:
    fq_prof_pacbio_ccs95 = default_pacbio_ccs95()
    fq_prof_pacbio_subreads = default_pacbio_subreads()
    fq_prof_illumina = default_illumina()

  ref = read_fasta_into_hash(args.reference_genome)
  txn = Transcriptome()
  txn.set_reference_genome_dictionary(ref)
  with open(args.transcripts_genepred) as inf:
    for line in inf:
      if line[0]=='#': continue
      txn.add_genepred_line(line.rstrip())
  if args.isoform_expression:
    sys.stderr.write("Reading expression from a TSV\n")
    with open(args.isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        txn.add_expression(f[0],float(f[1]))
  elif args.uniform_expression:
    sys.stderr.write("Using uniform expression model\n")
  elif args.cufflinks_isoform_expression:
    sys.stderr.write("Using cufflinks expression\n")
    with open(args.cufflinks_isoform_expression) as inf:
      line1 = inf.readline()
      for line in inf:
        f = line.rstrip().split("\t")
        txn.add_expression(f[0],float(f[9]))
  sys.stderr.write("have transcriptome\n")
  for n in txn.ref_hash.keys(): del txn.ref_hash[n]
  rbe = SimulationBasics.RandomTranscriptomeEmitter(txn)
  # Now we have the transcriptomes set
  #Now our dataset is set up
  if args.short_reads_only:
    rbe.set_gaussian_fragmentation_default_hiseq()
    for zi in range(0,args.short_read_count):
      [name,seq] = rbe.emit_short_read(args.short_read_length)
      if args.no_errors:
        print "@SRSIM"+str(zi+1)
        print seq
        print "+"
        print 'I'*len(seq)
      else:
        l1perm = fq_prof_illumina.create_fastq_and_permute_sequence(seq)
        print "@SRSIM"+str(zi+1)
        print l1perm['seq']
        print "+"
        print l1perm['qual']
    return
  if args.long_reads_only:
    rbe.set_gaussian_fragmentation_default_pacbio()
    for zi in range(0,args.long_read_count):
      [name,seq] = rbe.emit_long_read()
      if args.no_errors:
        g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(zi+1)+'/ccs'
        print "@"+g
        print seq
        print "+"
        print 'I'*len(seq)   
      else: 
        g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(zi+1)+'/ccs'
        seqperm = fq_prof_pacbio_ccs95.create_fastq_and_permute_sequence(seq)
        print "@"+g
        print seqperm['seq']
        print "+"
        print seqperm['qual']  
    return
  if not os.path.exists(args.output):
    os.makedirs(args.output)


  rbe.set_gaussian_fragmentation_default_hiseq()
  # Lets prepare to output now
  sys.stderr.write("Sequencing short reads\n")
  global left_handle
  global right_handle
  left_handle = gzip.open(args.output+"/SR_1.fq.gz",'wb')
  right_handle = gzip.open(args.output+"/SR_2.fq.gz",'wb')
  buffer_size = 10000
  buffer = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  z = 0
  for i in range(0,args.short_read_count):
    z = i+1
    if z %1000==0: sys.stderr.write(str(z)+"\r")
    buffer.append(z)
    if len(buffer) >= buffer_size:
      if args.threads <= 1:
        v = process_short_read_buffer(buffer[:],rbe,args,fq_prof_illumina)
        do_short(v)
      else:
        p.apply_async(process_short_read_buffer,args=(buffer[:],rbe,args,fq_prof_illumina),callback=do_short)
      buffer = []
  if len(buffer) > 0:
    if args.threads <= 1:
      v = process_short_read_buffer(buffer[:],rbe,args,fq_prof_illumina)
      do_short(v)
    else:
      p.apply_async(process_short_read_buffer,args=(buffer[:],rbe,args,fq_prof_illumina),callback=do_short)
    buffer = []
  if args.threads > 1:
    p.close()
    p.join()

  global greport
  of = open(args.output+"/SR_report.txt",'w')
  for name in greport:
    of.write("\t".join([str(x) for x in greport[name]])+"\n")
  of.close()  
  greport = {}

  sys.stderr.write("\nFinished sequencing short reads\n")
  left_handle.close()
  right_handle.close()

  # Now lets create the long read set
  rbe.set_gaussian_fragmentation_default_pacbio()
  sys.stderr.write("Sequencing ccs long reads\n")
  global long_handle
  long_handle = gzip.open(args.output+"/LR_ccs.fq.gz",'wb')
  buffer_size = 1000
  buffer = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  z = 0
  for i in range(0,args.long_read_count):
    z = i+1
    if z %100==0: sys.stderr.write(str(z)+"\r")
    buffer.append(z)
    if len(buffer) >= buffer_size:
      if args.threads <= 1:
        v = process_long_reads(buffer[:],rbe,args,fq_prof_pacbio_ccs95,'ccs')
        do_long(v)
      else:
        p.apply_async(process_long_reads,args=(buffer[:],rbe,args,fq_prof_pacbio_ccs95,'ccs'),callback=do_long)
      buffer = []
  if len(buffer) > 0:
    if args.threads <= 1:
      v = process_long_reads(buffer[:],rbe,args,fq_prof_pacbio_ccs95,'ccs')
      do_long(v)
    else:
      p.apply_async(process_long_reads,args=(buffer[:],rbe,args,fq_prof_pacbio_ccs95,'ccs'),callback=do_long)
    buffer = []
  if args.threads > 1:
    p.close()
    p.join()

  long_handle.close()
  of = open(args.output+"/LR_ccs_report.txt",'w')
  for name in greport:
    of.write("\t".join([str(x) for x in greport[name]])+"\n")
  of.close()  
  greport = {}
  sys.stderr.write("\nFinished sequencing ccs long reads\n")

  sys.stderr.write("Sequencing long sub reads\n")
  long_handle = gzip.open(args.output+"/LR_sub.fq.gz",'wb')
  buffer_size = 1000
  buffer = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for i in range(z,z+args.long_read_count):
    z = i+1
    if z %100==0: sys.stderr.write(str(z)+"\r")
    buffer.append(z)
    if len(buffer) >= buffer_size:
      if args.threads <= 1:
        v = process_long_reads(buffer[:],rbe,args,fq_prof_pacbio_subreads,'sub')
        do_long(v)
      else:
        p.apply_async(process_long_reads,args=(buffer[:],rbe,args,fq_prof_pacbio_subreads,'sub'),callback=do_long)
      buffer = []
  if len(buffer) > 0:
    if args.threads <= 1:
      v = process_long_reads(buffer[:],rbe,args,fq_prof_pacbio_subreads,'sub')
      do_long(v)
    else:
      p.apply_async(process_long_reads,args=(buffer[:],rbe,args,fq_prof_pacbio_subreads,'sub'),callback=do_long)
    buffer = []
  if args.threads > 1:
    p.close()
    p.join()

  long_handle.close()
  of = open(args.output+"/LR_sub_report.txt",'w')
  for name in greport:
    of.write("\t".join([str(x) for x in greport[name]])+"\n")
  of.close()  
  greport = {}
  sys.stderr.write("\nFinished sequencing long sub reads\n")

  combo = {}
  with open(args.output+"/SR_report.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      [name,express,left] = f
      if name not in combo:
        combo[name] = {}
        combo[name]['express'] = express
        combo[name]['left'] = 0
      combo[name]['left'] += int(left)
  with open(args.output+"/LR_ccs_report.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      [name,express,left] = f
      if name not in combo:
        combo[name] = {}
        combo[name]['express'] = express
        combo[name]['left'] = 0
      combo[name]['left'] += int(left)
  with open(args.output+"/LR_sub_report.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      [name,express,left] = f
      if name not in combo:
        combo[name] = {}
        combo[name]['express'] = express
        combo[name]['left'] = 0
      combo[name]['left'] += int(left)
  of = open(args.output+"/LR_SR_combo_report.txt",'w')
  for name in sorted(combo):
    of.write(name+"\t"+combo[name]['express']+"\t"+str(combo[name]['left'])+"\n")
  of.close()

  
def process_long_reads(buffer,rbe,args,error_profile,name_type):
  long = ''
  for z in buffer:
    [name,seq] = rbe.emit_long_read()
    g = ''
    if name_type == 'ccs':
      g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/ccs'
    elif name_type == 'sub':
      g = 'm150101_010101_11111_c111111111111111111_s1_p0/'+str(z)+'/0_'+str(len(seq)-1)
    else:
      sys.stderr.write("unknown name type "+name_type+"\n")
      sys.exit()
    if args.no_errors:
      long += "@"+g+"\n"
      long += seq+"\n"
      long += "+\n"
      long += len(seq)*'I'+"\n"
    else:
      seqperm = error_profile.create_fastq_and_permute_sequence(seq)
      long += "@"+g+"\n"
      long += seqperm['seq']+"\n"
      long += "+\n"
      long += seqperm['qual']+"\n"
  report = []
  # Now lets print out some of the emission details
  for name in sorted(rbe.emissions_report.keys()):
    express = 1
    if rbe.transcriptome.expression:
      express = rbe.transcriptome.expression.get_expression(name)
    entry = [name,express,rbe.emissions_report[name]]
    report.append(entry)
  rbe.emissions_report = {}
  return [report,long]

def do_long(v):
  global long_handle
  global glock
  if not v: return
  [report, long] = v
  glock.acquire()
  long_handle.write(long)
  add_report(report)
  glock.release()  
  return


def process_short_read_buffer(buffer,rbe,args,fq_prof_illumina):
  left = ''
  right = ''
  for z in buffer:
    [name,l1,r1] = rbe.emit_paired_short_read(args.short_read_length)
    if args.no_errors:
      left+="@SRSIM"+str(z)+"\n"
      left+=l1+"\n"
      left+="+\n"
      left+=len(l1)*'I'+"\n"
      right+="@SRSIM"+str(z)+"\n"
      right+=r1+"\n"
      right+="+\n"
      right+=len(r1)*'I'+"\n"
    else:
      left+="@SRSIM"+str(z)+"\n"
      l1perm = fq_prof_illumina.create_fastq_and_permute_sequence(l1)
      left+=l1perm['seq']+"\n"
      left+="+\n"
      left+=l1perm['qual']+"\n"
      r1perm = fq_prof_illumina.create_fastq_and_permute_sequence(r1)
      right+="@SRSIM"+str(z)+"\n"
      right+=r1perm['seq']+"\n"
      right+="+\n"
      right+=r1perm['qual']+"\n"

  report = []
  ## Now lets print out some of the emission details
  #of = open(args.output+"/SR_report.txt",'w')
  for name in sorted(rbe.emissions_report.keys()):
    express = 1
    if rbe.transcriptome.expression:
      express = rbe.transcriptome.expression.get_expression(name)
    entry = [name,express,rbe.emissions_report[name]]
    report.append(entry)
    #of.write(name +"\t"+str(express)+"\t"+str(rbe.emissions_report[name])+"\n")
  rbe.emissions_report = {}
  return [report,left,right]

def do_short(v):
  global left_handle
  global right_handle
  global glock
  if not v: return
  [report, left, right] = v
  glock.acquire()
  left_handle.write(left)
  right_handle.write(right)
  add_report(report)
  glock.release()  

def add_report(report):
  global greport
  for e in report:
    if e[0] not in greport:
      greport[e[0]] = e
    else:
      greport[e[0]][2] += e[2]
  return

if __name__=="__main__":
  main()
