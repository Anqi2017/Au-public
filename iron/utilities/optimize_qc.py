#!/usr/bin/python
import sys, os, subprocess, multiprocessing, re, zlib, argparse
import SamBasics
from SequenceBasics import GenericFastqFileReader, read_fasta_into_hash
from random import randint
from shutil import rmtree

# Test qc decisions on a fastq file for parameters such as
# (1) Left trim
# (2) Right trim
# (3) Minimum tolerated quality score for (4) X number of bases
# (5) Maximum tolerated number of mismatches in mapped bases

def main():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  cpus = multiprocessing.cpu_count()
  parser.add_argument('--ref_genome',required=True,help='(required) FASTA filename of reference genome')
  parser.add_argument('--bwa_index',required=True,help='(required) BWA Index')
  #parser.add_argument('--max_mismatches',type=int,default=2,help='INT maximum number of allowed mismatches')
  parser.add_argument('--min_read_size',type=int,default=30,help='INT minimum read size to consider')
  parser.add_argument('--test_size',type=int,default=5000,help='INT number of sequences to test')
  parser.add_argument('--min_test_size',type=int,default=500,help='INT disregard any parameter sets that do not produce at least this number of sequences prior to mapping.')
  parser.add_argument('--left_trim_range',help='start:end:increment, default is 0:[read_length]:5')
  parser.add_argument('--right_trim_range',help='start:end:increment, default is 0:[read_length]:5')
  parser.add_argument('--quality_number_range',help='start:end:increment, default is [qual_min]:[qual_max]:5')
  parser.add_argument('--quality_fail_count_range',help='start:end:increment, default is 0:[read_length]:5')
  parser.add_argument('--mapped_mismatch_range',help='start:end:increment, default is 0:3:1')
  parser.add_argument('--ignore_mapped_mismatches',action='store_true')
  parser.add_argument('--ignore_quality',action='store_true')
  parser.add_argument('--threads',type=int,default=cpus,help='INT of threads to run that defaults to cpu_count')
  parser.add_argument('--tempdir',default='/tmp/',help='Directory of your prefered temporary directory')
  parser.add_argument('-o',help='FILENAME for output')
  parser.add_argument('fastq_file',help='FILENAME for fastq file (can be .gz)')
  args = parser.parse_args()
  maxcnt = args.test_size
  mincnt = args.min_test_size
  sys.stderr.write("Testing up to "+str(maxcnt)+" reads.\n")
  sys.stderr.write("Require parameters leave at least "+str(mincnt)+" reads.\n")
  #max_allowed_mismatches = args.max_mismatches
  #sys.stderr.write("Allowing up to "+str(max_allowed_mismatches)+" mismatches.\n")
  #max_end_mismatches = 2
  min_read_size = args.min_read_size
  sys.stderr.write("Requiring QC parameters produce a minimum read length of "+str(min_read_size)+"\n")
  man = multiprocessing.Manager()
  Q = man.Queue()
  ifile = args.bwa_index 
  sys.stderr.write("BWA index: "+ifile+"\n")
  refgenome = args.ref_genome
  sys.stderr.write("Ref Genome: "+refgenome+"\n")
  #ifile = '/Shared/Au/jason/Reference/UCSC/Human/hg19_GRCh37_feb2009/BWA_Index/genome.fa'
  #refgenome = '/Users/weirathe/jason/Reference/UCSC/Human/hg19_GRCh37_feb2009/Genome/genome.fa'
  #refgenome = 'test_ref.fa'
  if args.threads: cpus = args.threads
  sys.stderr.write("Using "+str(cpus)+" threads\n")
  sys.stderr.write("reading reference genome\n")
  g = read_fasta_into_hash(refgenome)
  gz = {}
  cman = multiprocessing.Manager()
  cQ = man.Queue()
  pc = multiprocessing.Pool(processes=cpus)
  cresults = []
  sys.stderr.write("compressing reference genome\n")
  for name in g:
    pc.apply_async(comp,[name,g[name],cQ,len(g)])
  pc.close()
  pc.join()
  sys.stderr.write("\n")
  while not cQ.empty():
    [name,zseq] = cQ.get()
    gz[name] = zseq

  sys.stderr.write("finished processing reference genome\n")

  #[entries,stats] = read_fastq('test3.fq',maxcnt)
  [entries,stats] = read_fastq(args.fastq_file,maxcnt)

  #tstart = '/tmp'
  tstart = args.tempdir.rstrip('/')
  tdir = tstart.rstrip('/')+'/'+'weirathe.'+str(randint(1,100000000))
  if not os.path.exists(tdir): os.makedirs(tdir)
  z = 0
  #max_l_trim = 20
  #max_r_trim = 20
  max_l_trim = stats['lenmax']
  max_r_trim = stats['lenmax']
  min_l_trim = 0
  min_r_trim = 0
  l_trim_iter = 5
  r_trim_iter = 5
  if args.left_trim_range:
    m = re.match('(\d+):(\d+):(\d+)',args.left_trim_range)
    if not m:
      sys.stderr.write("Error. malformed left trim range "+args.left_trim_range+"\n")
      return
    max_l_trim = int(m.group(2))
    min_l_trim = int(m.group(1))
    l_trim_iter = int(m.group(3))
  if args.right_trim_range:
    m = re.match('(\d+):(\d+):(\d+)',args.right_trim_range)
    if not m:
      sys.stderr.write("Error. malformed right trim range "+args.right_trim_range+"\n")
      return
    max_r_trim = int(m.group(2))
    min_r_trim = int(m.group(1))
    r_trim_iter = int(m.group(3))

  max_q_num = stats['qmax']
  max_q_fail = stats['lenmax']
  min_q_num = stats['qmin']
  min_q_fail = 0
  q_num_iter = 5
  q_fail_iter = 5

  if args.quality_number_range:
    m = re.match('(\d+):(\d+):(\d+)',args.quality_number_range)
    if not m:
      sys.stderr.write("Error. malformed quality number range "+args.quality_number_range+"\n")
      return
    max_q_num = int(m.group(2))
    min_q_num = int(m.group(1))
    q_num_iter = int(m.group(3))
  if args.quality_fail_count_range:
    m = re.match('(\d+):(\d+):(\d+)',args.quality_fail_count_range)
    if not m:
      sys.stderr.write("Error. malformed quality number range "+args.quality_fail_count_range+"\n")
      return
    max_q_fail = int(m.group(2))
    min_q_fail = int(m.group(1))
    q_fail_iter = int(m.group(3))

  if args.ignore_quality:
    max_q_fail = stats['lenmax']
    min_q_fail = stats['lenmax']
    q_fail_iter = 1
    max_q_num = stats['qmax']
    min_q_num = stats['qmax']
    q_num_iter = 1

  max_mismatch = 3
  min_mismatch = 0
  mismatch_iter = 1

  if args.mapped_mismatch_range:
    m = re.match('(\d+):(\d+):(\d+)',args.mapped_mismatch_range)
    if not m:
      sys.stderr.write("Error. malformed mapped mismatch tolerance range "+args.mapped_mismatch_range+"\n")
      return
    max_mismatch = int(m.group(2))
    min_mismatch = int(m.group(1))
    q_mismatch = int(m.group(3))

  if args.ignore_mapped_mismatches:
    min_mismatch = stats['lenmax']
    max_mismatch = stats['lenmax']
    mismatch_iter = 1

  flist = []
  run_params = {}
  run_stats = {}
  sys.stderr.write("Left trim search space: "+str(min_l_trim)+":"+str(min([stats['lenmax'],max_l_trim]))+":"+str(l_trim_iter)+"\n")
  sys.stderr.write("Right trim search space: "+str(min_r_trim)+":"+str(min([stats['lenmax'],max_r_trim]))+":"+str(r_trim_iter)+"\n")
  sys.stderr.write("Quality number search space: "+str(max(min_q_num,stats['qmin']))+":"+str(min(max_q_num,stats['qmax']))+":"+str(q_num_iter)+"\n")
  sys.stderr.write("Quality fail count search space: "+str(min_q_fail)+":"+str(min(stats['lenmax'],max_q_fail))+":"+str(q_fail_iter)+"\n")
  sys.stderr.write("Max mapped mismatch search space: "+str(min_mismatch)+":"+str(min(stats['lenmax'],max_mismatch))+":"+str(mismatch_iter)+"\n")
  for l_cut in range(min_l_trim,min([stats['lenmax'],max_l_trim])+1,l_trim_iter):
   for r_cut in range(min_r_trim,min([stats['lenmax'],max_r_trim])+1,r_trim_iter):
    for q_floor in range(max(min_q_num,stats['qmin']),min(max_q_num,stats['qmax'])+1,q_num_iter):
     for failure_limit in range(min(min_q_fail,stats['lenmax']-l_cut-r_cut),min(stats['lenmax']-l_cut-r_cut,max_q_fail)+1,q_fail_iter):
      for max_allowed_mismatches in range(min_mismatch,max_mismatch+1,mismatch_iter):
       z += 1
       run_params[z] = {}
       run_params[z]['l_cut'] = l_cut
       run_params[z]['r_cut'] = r_cut
       run_params[z]['q_floor'] = q_floor
       run_params[z]['failure_limit'] = failure_limit
       run_params[z]['max_allowed_mismatches'] = max_allowed_mismatches
       run_stats[z] = {}
       run_stats[z]['after_qc_reads'] = 0
       run_stats[z]['after_qc_bases'] = 0
       of = open(tdir+'/'+str(z)+'.fq','w')
       k = 0
       scnt = 0
       for e in entries:
         seq = e['seq']
         seq = left_trim(seq,l_cut)
         seq = right_trim(seq,r_cut)
         qual = e['quality']
         qual = left_trim(qual,l_cut)
         qual = right_trim(qual,r_cut)
         if len(seq) < min_read_size: continue
         failure_count = 0
         for i in range(0,len(qual)):
           if seq[i].upper() == 'N': failure_count += 1
           elif ord(qual[i]) < q_floor: failure_count += 1
         if failure_count > failure_limit: continue
         k+=1
         scnt += 1
         run_stats[z]['after_qc_reads'] += 1
         run_stats[z]['after_qc_bases'] += len(seq)
         of.write("@s_"+str(k)+"\n")
         of.write(seq+"\n")
         of.write('+'+"\n")
         of.write(qual+"\n")
       of.close()
       if scnt < mincnt: #how many sequences were left after filtering, make sure we have enough to care
         os.remove(tdir+'/'+str(z)+'.fq')
       else:
         flist.append(z)
  sys.stderr.write("total of "+str(len(flist))+" params\n")
  p = multiprocessing.Pool(processes=cpus)
  results = []
  for z in flist:
    p.apply_async(check_parameters,(z,gz,ifile,tdir,run_params[z]['max_allowed_mismatches'],Q,len(flist)))
    #check_parameters(z,gz,ifile,tdir,max_end_mismatches,max_allowed_mismatches,Q)
    #print str(map_bases) + "\t" + str(map_reads)
  p.close()
  p.join()
  sys.stderr.write("\n")
  run_results = {}
  while True:
    if Q.empty(): break
    [z, reads, bases] = Q.get()
    #[z, reads, bases] = result
    run_results[z] = {}
    run_results[z]['after_mapped_reads'] = reads    
    run_results[z]['after_mapped_bases'] = bases    

  header = "left_cut_count\tright_cut_count\tmin_quality_value\tmax_quality_failure_count\tmax_mapped_mismatch_count\toriginal_read_count\toriginal_base_count\tpost_qc_read_count\tpost_qc_base_count\tmapped_reads\tmapped_bases"
  if args.o:
    of = open(args.o,'w')
    of.write(header+"\n")
  else:
    print header
  for z in sorted(run_results.keys()):
    ostring =  str(run_params[z]['l_cut']) + "\t" + str(run_params[z]['r_cut']) + "\t" + \
               str(run_params[z]['q_floor']) + "\t" + str(run_params[z]['failure_limit']) + "\t"
    ostring += str(run_params[z]['max_allowed_mismatches']) + "\t"
    ostring += str(stats['readcount']) + "\t" + str(stats['basecount']) + "\t"
    ostring += str(run_stats[z]['after_qc_reads']) + "\t" + str(run_stats[z]['after_qc_bases']) + "\t" 
    ostring += str(run_results[z]['after_mapped_reads']) + "\t" + str(run_results[z]['after_mapped_bases']) + "\t" 
    if args.o:
      of.write(ostring+"\n")
    else:
      print ostring
  if args.o:
    of.close()
  rmtree(tdir)

def comp(name,seq,cQ,tot):
  res = [name, zlib.compress(seq.upper())]
  cQ.put(res)
  sys.stderr.write('\r'+(' '*30))
  sys.stderr.write('\r'+str(cQ.qsize())+'/'+str(tot))
  sys.stderr.flush()
  return 
  #[name, zlib.compress(seq.upper())]

def check_parameters(z,gz,ifile,tdir,max_allowed_mismatches,Q,fsize):
    #sys.stderr.write("doing "+str(z)+"\n")
    g = {}
    for n in gz:
      g[n] = zlib.decompress(gz[n])
    FNULL = open(os.devnull,'w')
    cmd1 = "bwa mem "+ifile+" "+tdir+'/'+str(z)+'.fq'
    cmd2 = "samtools view -S -"
    stream1 = subprocess.Popen(cmd1.split(),stdout=subprocess.PIPE,stderr=FNULL)
    stream2 = subprocess.Popen(cmd2.split(),stdin=stream1.stdout,stdout=subprocess.PIPE,stderr=FNULL)
    reads = {}
    while True:
      sumlen= 0
      mismatches = 0
      line = stream2.stdout.readline()
      if not line: break
      f = line.rstrip().split("\t")
      if f[2] == '*':
        continue
      d = SamBasics.sam_line_to_dictionary(line)
      #if d['rname'] != 'chr20': continue #get rid of this line soon.
      cigar = d['cigar_array']
      #endmismatch = 0
      #if cigar[0]['op'] == 'S':
      #  endmismatch += cigar[0]['val']
      #if cigar[len(cigar)-1]['op'] == 'S':
      #  endmismatch += cigar[len(cigar)-1]['val']
      #if endmismatch > max_end_mismatches: continue
      read_index = 1
      chrom_index = d['pos']
      for e in cigar:
        if re.match('[MX=]',e['op']): 
          sumlen += e['val']  # keep track of our match length
          refseq = g[d['rname']][chrom_index-1:chrom_index-1+e['val']].upper()
          readseq = d['seq'][read_index-1:read_index-1+e['val']].upper()
          for i in range(0,e['val']): 
            if refseq[i] != readseq[i]: mismatches += 1
          read_index += e['val']
          chrom_index += e['val']
        elif re.match('[SI]',e['op']):
          mismatches += e['val']
          read_index += e['val']
        elif re.match('[NDH]',e['op']):
          chrom_index += e['val']
        else:
          sys.stderr.write("warning: strange SAM op\n")
      # save the biggest sum for the read name
      #print 'mismatches: '+str(mismatches)
      if mismatches > max_allowed_mismatches: continue
      if d['qname'] not in reads: 
        reads[d['qname']] = {}
        reads[d['qname']]['alignment_length'] = 0
        reads[d['qname']]['mismatches'] = 0
      if sumlen > reads[d['qname']]['alignment_length']: 
        reads[d['qname']]['alignment_length'] = sumlen
        reads[d['qname']]['mismatches'] = mismatches
    mapped_bases = 0
    mapped_reads = 0
    for rname in reads:
      mapped_bases += reads[rname]['alignment_length']
      mapped_reads += 1
    #print str(mapped_bases) + "\t" + str(mapped_reads)
    res = [z,mapped_reads,mapped_bases]
    #sys.stderr.write(str(z)+"\t"+str(mapped_reads)+"\t"+str(mapped_bases)+"\n")
    Q.put(res)
    progress = Q.qsize()
    sys.stderr.write('\r'+(' '*40))
    sys.stderr.write('\r'+str(progress)+"/"+str(fsize))
    sys.stderr.flush()
    return

def read_fastq(fastq_file,maxcnt):
  gfr = GenericFastqFileReader(fastq_file)
  ecnt = 0
  qseen = set()
  lenmax = 0
  lenmin = float('inf')
  entries = []
  bases = 0
  while True:
    e = gfr.read_entry()
    if not e or ecnt > maxcnt: break
    ecnt += 1
    slen = len(e['seq'])
    if slen < lenmin: lenmin = slen
    if slen > lenmax: lenmax = slen
    seq = e['seq']
    bases += len(seq)
    for v in [ord(x) for x in e['quality']]:
      qseen.add(v)
    entries.append(e)
  gfr.close()
  qmin = min(qseen)
  qmax = max(qseen)
  stats  = {}
  stats['qmin'] = qmin
  stats['qmax'] = qmax
  stats['lenmin'] = lenmin
  stats['lenmax'] = lenmax
  stats['readcount'] = len(entries)
  stats['basecount'] = bases
  return [entries,stats]

def right_trim(seq,n):
  if n == 0: return seq
  return seq[:-n]

def left_trim(seq,n):
  if n == 0: return seq
  return seq[n:]  

main()
