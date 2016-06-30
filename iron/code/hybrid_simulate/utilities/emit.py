#!/usr/bin/python
import argparse, sys, os, pickle, zlib, base64, json, math, gzip, re
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from subprocess import PIPE, Popen
from Bio.Simulation.Emitter import TranscriptomeEmitter
from Bio.Structure import Transcriptome
from Bio.Simulation.RandomSource import RandomSource
from Bio.Simulation.Permute import MakeCuts, random_flip, MakeErrors, rate_to_phred33, phred33_to_rate
from Bio.Sequence import rc
from Bio.Format.Fastq import Fastq

def main(args):
  # check outputs
  if len(args.output) > 1 and not args.sr:
    sys.stderr.write("Error: Long reads don't support multiple output files\n")
    sys.exit()
  elif len(args.output) > 2:
    sys.stderr.wrtie("Error: Short reads support at most two output files (paired end)\n")
    sys.exit()
  if args.sr_length < args.minimum_read_length:
    args.minimum_read_length = args.sr_length
  inf = sys.stdin
  if args.emitter != '-':
    inf = open(args.emitter)
  sys.stderr.write("reading in transcriptome emitter\n")
  indata = pickle.loads(zlib.decompress(base64.b64decode(inf.read().rstrip())))
  txome = Transcriptome()
  txome.load_serialized(indata['txome'])
  rnum = RandomSource()
  rnum_tx = RandomSource() # for drawing transcripts
  if args.seed: 
    rnum = RandomSource(args.seed)
    rnum_tx = RandomSource(args.seed)
  elif indata['seed']: 
    rnum = RandomSource(args.seed)
    rnum_tx = RandomSource(args.seed)
  # Load in error profile data
  ep = None
  if args.error_profile:
    sys.stderr.write("read in error profile\n")
    ep = ErrorProfilePermuter(args.error_profile,rnum,args.skew_profile_error_rate)
  txemitter = TranscriptomeEmitter(txome,rand=rnum_tx)
  if indata['weight_type'] == 'expression_table':
    sys.stderr.write("Using expression table defined transcript expression\n")
    txweight = indata['weights']
    txemitter.set_weights_by_dict(txweight)
  elif indata['weight_type'] == 'exponential_distribution':
    sys.stderr.write("ERROR not yet implemented exponential distribution\n")
    sys.exit()
  elif indata['weight_type'] == 'uniform_distribution':
    sys.stderr.write("Using uniform distribution of transcript expression\n")
  cutter = MakeCuts(rand=rnum_tx)
  if args.sr:
    cutter.set_custom(args.sr_gauss_min,args.sr_gauss_mu,args.sr_gauss_sigma)
  elif args.lr:
    cutter.set_custom(args.lr_gauss_min,args.lr_gauss_mu,args.lr_gauss_sigma)
  # Prepare outputs
  of1 = sys.stdout
  if args.output[0][-3:] == '.gz':
    of1 = gzip.open(args.output[0],'w')
  elif args.output[0] != '-':
    of1 = open(args.output[0],'w')
  of2 = None
  if len(args.output) > 1:
    if args.output[1][-3:] == '.gz':
      of2 = gzip.open(args.output[1],'w')
    elif args.output[0] != '-':
      of2 = open(args.ouptput[1],'w')
  of_origin = None
  if args.output_original_source:
    if args.output_original_source[-3:]=='.gz':
      of_origin = gzip.open(args.output_original_source,'w')
    else:
      of_origin = open(args.output_original_source,'w')
  of_sc = None
  if args.output_sequence_change:
    if args.output_sequence_change[-3:]=='.gz':
      of_sc = gzip.open(args.output_sequence_change,'w')
    else:
      of_sc = open(args.output_sequence_change,'w')
  
  absmax = args.count*100
  finished_count = 0
  z = 0
  while finished_count < args.count:
    z += 1
    if z > absmax: break
    tx = txemitter.emit_transcript()
    seq = tx.get_sequence()
    stage1seq = seq
    if args.trim_5prime or args.trim_3prime:
      fivestart = 0
      threeend = len(seq)
      if args.trim_5prime:
        lcut = int(args.trim_5prime[0])*len(seq)
        rcut = int(args.trim_5prime[1])*len(seq)
        fivestart = rnum_tx.randint(lcut,rcut)
      if args.trim_3prime:
        lcut = int(args.trim_3prime[0])*len(seq)
        rcut = int(args.trim_3prime[1])*len(seq)
        threeend = rnum_tx.randint(lcut,rcut)
      # set sequence to its new trimmed bounds
      seq = seq[fivestart:threeend]

    # flip sequence if necessary
    if not args.no_flip:
      seq = random_flip(seq,rnum_tx)

    l_read = create_name(rnum)
    r_read = None
    if args.sr or args.lr:
     cutseq = cutter.get_cut(seq)
    else: cutseq = seq #case for no_fragmentation
    ############# if we pass this we will really start with this one
    if len(cutseq) < args.minimum_read_length: continue
    # can now log our read name
    if of_origin:
      of_origin.write(l_read+"\t"+tx.get_gene_name()+"\t"+tx.get_transcript_name()+"\n")
    stage2seq = cutseq
    r = None
    if args.sr:
      r_read = l_read
      l = cutseq[0:args.sr_length]
      r = rc(cutseq[-1*args.sr_length:])
    elif args.lr:
      l = cutseq
    else: l = cutseq
    stage3left = l
    stage3right = r
    if not stage3right: stage3right = ''
    #################
    #  l (or l and r) contains the sequence prior to errors being added
    l_qual = 'I'*len(l) 
    r_qual = None
    if r: r_qual = 'I'*len(r)
    if args.fixed_quality:
      sys.stderr.write("Use fixed quality\n")
      if len(args.fixed_quality) != 1:
        sys.stderr.write("ERROR fixed quaility should be 1 character\n")
        sys.exit()
      l_qual = args.fixed_quality*len(l)
      if r: r_qual = args.fixed_quality*len(r)
    elif args.quality_from_error_rate:
      #sys.stderr.write("Set quality from error rate\n")
      qchar = chr(int(-10*math.log10(args.quality_from_error_rate))+33)
      l_qual = qchar*len(l)
      if r: r_qual = qchar*len(r)
    else: #default is generate quality from profile
      if not ep:
        sys.stderr.write("ERROR: cannot generate quality from a profile.  Set error profile or chooce quaility from error rate or fixed quality\n")
        sys.exit()
      l_qual = ep.emit_qual(len(l))
      if r: r_qual = ep.emit_qual(len(r))
    # Now prior to errors l_qual and r_qual contain our qualities

    l_fastq = Fastq([l_read,l,'+',l_qual])
    r_fastq = None
    if r:
      r_fastq = Fastq([r_read,r,'+',r_qual])
    # Permute sequences by a specific error rate
    if args.specific_errors:
      rate = args.specific_errors
      me = MakeErrors(rand=rnum)
      if args.specific_before_context: me.set_before_context(args.specific_before_context)
      if args.specific_after_context: me.set_after_context(args.specific_after_context)
      if args.specific_reference_base: 
        if args.specific_reference_base != '-':
          me.set_observed_base(args.specific_reference_base)
      if args.specific_modified_base: 
        if args.specific_modified_base != '-':
          me.set_modified_base(args.specific_modified_base)
      if args.specific_reference_base == '-': #doing insertions
        l_fastq = me.random_insertion(l_fastq,rate)
        if r_fastq: r_fastq = me.random_insertion(r_fastq,rate)
      elif args.specific_modified_base == '-': #doing deletions
        l_fastq = me.random_deletion(l_fastq,rate)
        if r_fastq: r_fastq = me.random_insertion(r_fastq,rate)
      else:
        l_fastq = me.random_substitution(l_fastq,rate)
        if r_fastq: r_fastq = me.random_insertion(r_fastq,rate)
    elif args.uniform_any_error:
      l_fastq = do_uniform_any(l_fastq,rnum,args.uniform_any_error)
      if r_fastq: r_fastq = do_uniform_any(r_fastq,rnum,args.uniform_any_error)  
    elif args.uniform_mismatch_error:
      l_fastq = do_uniform_mismatch(l_fastq,rnum,args.uniform_mismatch_error)
      if r_fastq: r_fastq = do_uniform_mismatch(r_fastq,rnum,args.uniform_mismatch_error)  
    elif args.any_error_by_quality:
      l_fastq = do_quality_any(l_fastq,rnum)
      if r_fastq: r_fastq = do_quality_any(r_fastq,rnum)      
    elif args.mismatch_error_by_quality:
      l_fastq = do_quality_mismatch(l_fastq,rnum)
      if r_fastq: r_fastq = do_quality_mismatch(r_fastq,rnum)
    elif args.profile_context_error:
      l_fastq = ep.permute_context(l_fastq)
      if r_fastq: r_fastq = ep.permute_context(r_fastq)
    elif args.profile_general_error:
      l_fastq = ep.permute_general(l_fastq)
      if r_fastq: r_fastq = ep.permute_general(r_fastq)
      
    # if SR grown/shrink to appropriate length
    if args.sr and len(l_fastq) != args.sr_length:
      l_fastq = fit_length(l_fastq,args.sr_length,rnum)
    if r:
      if args.sr and len(r_fastq) != args.sr_length:
        r_fastq = fit_length(r_fastq,args.sr_length,rnum)

    of1.write(l_fastq.fastq())
    if of2: 
      of2.write(r_fastq.fastq())

    stage4left = l_fastq.seq
    stage4right = ''
    if of_sc:
      of_sc.write(l_fastq.name+"\t"+tx.get_gene_name()+"\t"+tx.get_transcript_name()+"\t" \
                + stage1seq+"\t"+stage2seq+"\t"+stage3left+"\t"+stage3right+"\t"+stage4left+"\t"+stage4right+"\n")
    if r_fastq: stage4right = r_fastq.seq
    finished_count += 1
    sys.stderr.write(str(finished_count)+'/'+str(args.count)+"   \r")
  sys.stderr.write("\n")
  of1.close()
  if of2:
    of2.close()
  if of_origin:
    of_origin.close()
  if of_sc:
    of_sc.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_quality_mismatch(fastq,rnum):
  sequence = fastq.seq
  quality = fastq.qual
  #five possible errors 3 base changes a deletion or insertion
  seq = ''
  qual = ''
  for i in range(0,len(sequence)):
    rate = phred33_to_rate(quality[i])
    if rnum.random() > rate:
      seq += sequence[i]
      qual += quality[i]
      continue
    seq += rnum.different_random_nt(sequence[i])
    qual += quality[i]
    #for del type we do nothing ... cause its not getting added
  return Fastq([fastq.name,seq,'+',qual])

def do_quality_any(fastq,rnum):
  sequence = fastq.seq
  quality = fastq.qual
  #five possible errors 3 base changes a deletion or insertion
  seq = ''
  qual = ''
  for i in range(0,len(sequence)):
    rate = phred33_to_rate(quality[i])
    if rnum.random() > rate:
      seq += sequence[i]
      qual += quality[i]
      continue
    type = rnum.choice(['ins','del','mis','mis','mis'])
    if type == 'mis':
      seq += rnum.different_random_nt(sequence[i])
      qual += quality[i]
    elif type == 'ins':
      if rnum.random() < 0.5:
        seq += sequence[i]+rnum.random_nt()
        qual += quality[i]+quality[i]
      else:
        seq += rnum.random_nt()+sequence[i]
        qual += quality[i]+quality[i]
    #for del type we do nothing ... cause its not getting added
  return Fastq([fastq.name,seq,'+',qual])

def do_uniform_any(fastq,rnum,rate):
  #ibase = rate_to_phred33(rate)
  #six error possibilies
  sequence = fastq.seq
  quality = fastq.qual
  #five possible errors 3 base changes a deletion or insertion
  seq = ''
  qual = ''
  for i in range(0,len(sequence)):
    if rnum.random() > rate:
      seq += sequence[i]
      qual += quality[i]
      continue
    type = rnum.choice(['ins','del','mis','mis','mis'])
    if type == 'mis':
      seq += rnum.different_random_nt(sequence[i])
      qual += quality[i]
    elif type == 'ins':
      if rnum.random() < 0.5:
        seq += sequence[i]+rnum.random_nt()
        qual += quality[i]+quality[i]
      else:
        seq += rnum.random_nt()+sequence[i]
        qual += quality[i]+quality[i]
    #for del type we do nothing ... cause its not getting added
  return Fastq([fastq.name,seq,'+',qual])

def do_uniform_mismatch(fastq,rnum,rate):
  #ibase = rate_to_phred33(rate)
  #six error possibilies
  sequence = fastq.seq
  quality = fastq.qual
  #five possible errors 3 base changes a deletion or insertion
  seq = ''
  qual = ''
  for i in range(0,len(sequence)):
    if rnum.random() > rate:
      seq += sequence[i]
      qual += quality[i]
      continue
    seq += rnum.different_random_nt(sequence[i])
    qual += quality[i]
    #for del type we do nothing ... cause its not getting added
  return Fastq([fastq.name,seq,'+',qual])

def fit_length(fastq,target_length,rnum):
  sequence = fastq.seq
  quality = fastq.qual
  if len(sequence)==target_length: return fastq
  bps = ['A','C','G','T']
  side = ['left','right']
  ibase = rate_to_phred33(0.75)
  if len(sequence) < target_length:
    while len(sequence) < target_length:
      rside = rnum.choice(side)
      if rside == 'left':
        sequence = rnum.choice(bps)+sequence
        quality  = ibase+quality
      else:
        sequence = sequence+rnum.choice(bps)
        quality = quality+ibase
    return Fastq([fastq.name,sequence,'+',quality])
  # must be bigger
  while len(sequence) > target_length:
    rside = rnum.choice(side)
    if rside == 'left':
      sequence = sequence[1:]
      quality = quality[1:]
    else:
      sequence = sequence[:-1]
      quality = quality[:-1]
  return Fastq([fastq.name,sequence,'+',quality])

def create_name(r):
  return rhex(8,r)+'-'+rhex(4,r)+'-'+rhex(4,r)+'-'+rhex(4,r)+'-'+rhex(12,r)

def rhex(hlen,rnum):
  vals = 'abcdef0123456789'
  #out = ''
  return ''.join([rnum.choice(vals) for i in range(hlen)])

class ErrorProfilePermuter:
  def __init__(self,fname,rnum,skew=None):
    invals = json.loads(zlib.decompress(base64.b64decode(open(fname).read().rstrip())))
    self.random = rnum
    self.quality_counts = []
    #expand our quality counts
    sys.stderr.write("unpack qualities\n")
    for vals in invals['quality_counts']:
      self.quality_counts.append([])
      for part in [[chr(x[0])*x[1]]*x[2] for x in vals]:
        self.quality_counts[-1]+=part
    sys.stderr.write("finished unpacking\n")
    self.context_error = invals['context_error']
    self.alignment_error = invals['alignment_error']
    self.error_stats = invals['error_stats']
    # get error rate from error stats
    err_count = [int(y[1]) for y in [x.split("\t") for x in self.error_stats.rstrip().split("\n")] if y[0]=='ANY_ERROR'][0]
    base_count = [int(y[1]) for y in [x.split("\t") for x in self.error_stats.rstrip().split("\n")] if y[0]=='ALIGNMENT_BASES'][0]
    self.error_rate = float(err_count)/float(base_count)
    sys.stderr.write("error profile has error rate of "+str(self.error_rate)+"\n")
    self.md = {} #mismatches and deltions
    self.ins = {}
    for entry in self.context_error['data']:
      b = entry[0]
      a = entry[1]
      r = entry[2]
      m = entry[3]
      f = entry[4]
      if r == '-': # insertion
        if b not in self.ins: self.ins[b] = {}
        if a not in self.ins[b]: self.ins[b][a] = {}
        self.ins[b][a][m] = f
      else: # mismatches and deletions
        if b not in self.md: self.md[b] = {}
        if a not in self.md[b]: self.md[b][a] = {}
        if r not in self.md[b][a]: self.md[b][a][r] = {}
        self.md[b][a][r][m] = f
    #loaded in context
    self.gins = {}
    self.gmd = {}
    rtot = {}
    ginstot = self.alignment_error['data'][0][3] # total number of positions insertions could occur at should be similar to total checked
    # get totals for each reference context
    for e in self.alignment_error['data']:
      r = e[0]
      m = e[1]
      if r != '-':
        if r not in rtot: rtot[r] = 0
        rtot[r] += e[2]
    for e in self.alignment_error['data']:
      r = e[0]
      m = e[1]
      if r != '-':
        if r not in self.gmd: self.gmd[r] = {}
        if rtot[r] == 0:
          self.gmd[r][m] = 0
        else:
          self.gmd[r][m] = float(e[2])/float(rtot[r])
      else:
        if m != '-': # we aren't presently tracking them so we'll just calculate them from the others
          self.gins[m] = float(e[2])/ginstot
    self.gins['-'] = 1 - sum([self.gins[x] for x in self.gins.keys() if x != '-'])
    #loaded in general
    # if we want to nudge the error rate closer to some particular value while maintaining general pattern
    if skew:
      sys.stderr.write("try to go from "+str(self.error_rate)+" to "+str(skew)+"\n")
      if self.error_rate == 0:
        sys.stderr.write("ERROR: no errors to skew\n")
        sys.exit()
      factor = skew/self.error_rate
      for b in self.md:
        for a in self.md[b]:
          for r in self.md[b][a]:
            for m in self.md[b][a][r]:
              if r != m: self.md[b][a][r][m] = self.md[b][a][r][m]*factor
            #adjust the match
            self.md[b][a][r][r] = 1 - min(sum([self.md[b][a][r][x] for x in self.md[b][a][r] if x != r]),1)
      for b in self.ins:
        for a in self.ins[b]:
          for m in self.ins[b][a]:
            if m != '-': self.ins[b][a][m] = self.ins[b][a][m]*factor
          self.ins[b][a]['-'] = 1 - min(sum([self.ins[b][a][x] for x in self.ins.keys() if x != '-']),1)
      for r in self.gmd:
        for m in self.gmd[r]:
          if r != m: self.gmd[r][m] = self.gmd[r][m]*factor
        self.gmd[r][r] = 1 - min(sum([self.gmd[r][x] for x in self.gmd[r] if x != r]),1)
      for m in self.gins:
        if m != '-': self.gins[m] = self.gins[m]*factor
      self.gins['-'] = 1 - min(sum([self.gins[x] for x in self.gmd[r] if x!='-']),1)
  def draw_general_error(self,reference):
    md = self.gmd[reference]
    mods = sorted(md.keys())
    vals = [md[x] for x in mods]
    ind = self.random.get_weighted_random_index(vals)
    return mods[ind]

  def draw_general_insert(self):
    mods = sorted(self.gins.keys())
    vals = [self.gins[x] for x in mods]
    #print mods
    #print vals
    ind = self.random.get_weighted_random_index(vals)
    #print ind
    if mods[ind] == '-': return None
    return mods[ind]

  # based on error probability 
  def draw_context_error(self,before,after,reference):
    md = self.md[before][after][reference]
    #print md
    mods = sorted(md.keys())
    vals = [md[x] for x in mods]
    #print vals
    ind = self.random.get_weighted_random_index(vals)
    return mods[ind]
  def draw_context_insert(self,before,after):
    ins = self.ins[before][after]
    mods = sorted(ins.keys())
    vals = [ins[x] for x in mods]
    ind = self.random.get_weighted_random_index(vals)
    if mods[ind] == '-': return None
    return mods[ind]

  # modify a sequence by general error rates
  def permute_general(self,fastq):
    seq = ''
    qual = ''
    for i in range(0,len(fastq.seq)):
      # see about inserts
      bseq = ''
      bqual = ''
      if self.random.random() < 0.5:
        err = self.draw_general_insert()
        if err:
          bseq = err
          bqual = fastq.qual[i]
      aseq = ''
      aqual = ''
      if self.random.random() < 0.5:
        err = self.draw_general_insert()
        if err:
          aseq = err
          aqual = fastq.qual[i]
      err = self.draw_general_error(fastq.seq[i])
      if err != '-':
        seq += bseq+err+aseq
        qual += bqual+fastq.qual[i]+aqual
      else:
        seq += bseq+aseq
        qual += bqual+aqual
    return Fastq([fastq.name,seq,'+',qual])

  # modify sequence by context
  def permute_context(self,fastq):
    if len(fastq.seq) < 2: return fastq
    seq = fastq.seq[0]
    qual = fastq.qual[0]
    for i in range(1,len(fastq.seq)-1):
      insbefore = ''
      qualbefore = ''
      if self.random.random() < 0.5:
        insdrawn = self.draw_context_insert(fastq.seq[i-1],fastq.seq[i])
        if insdrawn: 
          insbefore = insdrawn
          qualbefore = fastq.qual[i]
      insafter = ''
      qualafter = ''
      if self.random.random() < 0.5:
        insdrawn = self.draw_context_insert(fastq.seq[i],fastq.seq[i+1])
        if insdrawn: 
          insafter = insdrawn
          qualafter = fastq.qual[i]
      err = self.draw_context_error(fastq.seq[i-1],fastq.seq[i],fastq.seq[i+1])
      if err != '-':
          seq += insbefore+err+insafter
          qual += qualbefore+fastq.qual[i]+qualafter
      else:
        seq += insbefore+insafter
        qual += qualbefore+qualafter
    return Fastq([fastq.name,seq+fastq.seq[-1],'+',qual+fastq.qual[-1]])

  def emit_qual(self,slen):
    full_len = ''
    while len(full_len) < slen:
      curr = len(full_len)
      bin = int(100*float(curr)/float(slen))
      val = self.random.choice(self.quality_counts[bin])
      full_len += val
    return full_len[0:slen]

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('emitter',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('count',type=int,help="Number of reads to output")
  parser.add_argument('--error_profile',help="Use a read profile")

  parser.add_argument('--minimum_read_length',type=int,default=200,help="Minimum read length (is over-ridden by sr_length if it is smaller)")
  parser.add_argument('--seed',type=int,help="Set a seed. If seed has been set in an emitter, then you set it here, this one replace the old one.")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  group5 = parser.add_argument_group(title="Output options")
  group5.add_argument('--output_original_source',help="Attribute each read to its original transcript and gene\n<read> <gene> <transcript>")
  group5.add_argument('--output_sequence_change',help="Output a table of how each read was perterbed by errors\n<read> <gene> <transcript> <transcript sequence> <post fragmentation/flipping> <left read pre-error> <right read pre-error> <left read post-error> <right read post-error>")
  group5.add_argument('-o','--output',nargs='+',required=True,help="OUTPUTFILE or STDOUT if not set")

  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--lr',action='store_true',help="Simulate long reads from transcripts")
  group.add_argument('--sr',action='store_true',help="Simulate short reads from transcripts")
  group.add_argument('--no_fragmentation',action='store_true',help="No fragmentation.")

  group1 = parser.add_argument_group(title="LR params",description="Long read parameters.")
  group1.add_argument('--lr_gauss_min',type=float,default=1000,help="Minimum value")
  group1.add_argument('--lr_gauss_mu',type=float,default=4000,help="Average value")
  group1.add_argument('--lr_gauss_sigma',type=float,default=500,help="Standard deviation")

  group2 = parser.add_argument_group(title="SR params",description="Long read parameters.")
  group2.add_argument('--sr_length',type=int,default=100,help="Short read length")
  group2.add_argument('--sr_gauss_min',type=float,default=150,help="Minimum value")
  group2.add_argument('--sr_gauss_mu',type=float,default=290,help="Average value")
  group2.add_argument('--sr_gauss_sigma',type=float,default=290,help="Standard deviation")

  group3 = parser.add_argument_group(title="Error params",description="Options for permuting sequence with errors")
  mgroup3 = group3.add_mutually_exclusive_group()
  mgroup3.add_argument('--uniform_any_error',type=float,help="use a constant error rate.")
  mgroup3.add_argument('--uniform_mismatch_error',type=float,help="use a constant error rate, mismatches only")
  mgroup3.add_argument('--profile_general_error',action='store_true',help="use the error profiles error rate")
  mgroup3.add_argument('--profile_context_error',action='store_true',help="use context specific error rates")
  mgroup3.add_argument('--any_error_by_quality',action='store_true',help="introduce any error type by quality")
  mgroup3.add_argument('--mismatch_error_by_quality',action='store_true',help="introduce mismatch_errors by quality")
  mgroup3.add_argument('--specific_errors',type=float,help="RATE at which to introduce very specific errors")
  group3.add_argument('--skew_profile_error_rate',type=float,help="Try to adjust the error rate to this number")

  group5 = parser.add_argument_group(title="Specific errors",description="inject very specific errors")
  group5.add_argument('--specific_before_context',help="Before base")
  group5.add_argument('--specific_after_context',help="After base")
  group5.add_argument('--specific_reference_base',help="Reference base")
  group5.add_argument('--specific_modified_base',help="Modified base")

  group4 = parser.add_argument_group(title="Quality params",description="Options for producing the quality string")
  mgroup4 = group4.add_mutually_exclusive_group()
  mgroup4.add_argument('--quality_from_profile',action='store_true',help="Use the profile for quality.  This is default.")
  mgroup4.add_argument('--quality_from_error_rate',type=float,help="Make a new quality string from the error rate. Phred33 format")
  mgroup4.add_argument('--fixed_quality',help="Quality string is this character")

  parser.add_argument('--trim_5prime',nargs=2,type=float,help="Keep 3' end. trim at random location between these fractions 0.0 is 5' end and 1.0 is 3' end")
  parser.add_argument('--trim_3prime',nargs=2,type=float,help="Keep 5' end. trim at random location between these fractions 0.0 is 5' end and 1.0 is 3' end")


  parser.add_argument('--no_flip',action='store_true',help="Don't randomly flip transcripts.  Keep strand specificity.")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()

  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

def external_cmd(cmd,version=None):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main()
