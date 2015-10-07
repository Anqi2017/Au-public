#!/usr/bin/python
import os, argparse, sys, re, random, multiprocessing, subprocess, json
import SamBasics, PSLBasics, GenePredBasics
from SequenceBasics import read_fasta_into_hash
from shutil import rmtree


def main():
  parser = argparse.ArgumentParser(description="Annotate output by a reference genepred")
  parser.add_argument('sam_file',help='SAM or BAM file')
  parser.add_argument('--reference_genome',help='FASTA file')
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default='/tmp',help="DIRECTORY in which we can make directories")
  group.add_argument('--specific_tempdir',help="DIRECTORY to create")
  parser.add_argument('--use_secondary_alignments',action='store_true',help="Use secondary alignemnts as well as primary.")
  parser.add_argument('--threads',type=int,default=0,help="Number of threads to use. Use max if not specified.")
  args = parser.parse_args()
  tdir = None
  if args.specific_tempdir:
    tdir = args.specific_tempdir.rstrip('/')
    #if os.path.exists(tdir):
    #  sys.stderr.write("ERROR specific temporary directory already exists.\n")
    #  sys.exit()
    if not os.path.exists(tdir):
      os.makedirs(tdir)
  else:
    tdir = args.tempdir.rstrip('/')+'/weirathe.'+str(random.randint(1,100000000))
    if not os.path.exists(tdir):
      os.makedirs(tdir)
  args.tempdir = tdir #replace tempdir
  sys.stderr.write("Working in: "+args.tempdir+"\n")

  if args.threads == 0: args.threads = multiprocessing.cpu_count()
  sys.stderr.write("Working on "+str(args.threads)+" threads\n")

  if args.threads > 1:
    p = multiprocessing.Pool(args.threads)
  for i in range(0,args.threads):
    if args.threads > 1:
      p.apply_async(make_exons,args=[args,i,args.threads])
    else:
      make_exons(args,i,args.threads)
  if args.threads > 1:
    p.close()
    p.join()

  #Now we can parse the bed file
  reads = {}
  for i in range(0,args.threads):
    with open(args.tempdir+'/bedpart.'+str(i)+'.bed') as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        t = {}
        t['chrom'] = f[0]
        t['tstart'] = int(f[1]) #bed format
        t['tend'] = int(f[2])
        t['strand'] = f[3]
        read = f[4]
        t['aligned_bases'] = int(f[5])
        t['indel_bases'] = int(f[6])
        t['alignment_type'] = f[7]
        t['qstart'] = int(f[8])
        if read not in reads:
          reads[read] = []
        reads[read].append(t)
  if args.threads > 1:
    p = multiprocessing.Pool(args.threads)
  for read in reads:
    locus_groups = get_locus_group(args,reads[read])
    #glines = get_genepred_lines(json.dumps(locus_groups),read)
    if args.threads > 1:
      p.apply_async(get_genepred_lines,args=(json.dumps(locus_groups),read),callback=docallback)
    else:
      glines = get_genepred_lines(json.dumps(locus_groups),read)
      docallback(glines)
    #for gline in glines:
    #  print gline.rstrip()
  if args.threads > 1:
    p.close()
    p.join()

  if not args.specific_tempdir:
    rmtree(tdir)

def docallback(glines):
  for line in glines:
    print line.rstrip()

def get_genepred_lines(jlocus_groups,read):
    locus_groups = json.loads(jlocus_groups)
    glines = []
    for entries in locus_groups:
      sentries = sorted(entries,key=getKey)
      strand = sentries[0]['strand']
      sp = subprocess.Popen('bedtools sort -i - | bedtools merge -i - | bedtools sort -i',shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
      for sentry in sentries:
        sp.stdin.write(sentry['chrom']+"\t"+str(sentry['tstart'])+"\t"+str(sentry['tend'])+"\n")
      mlines = [x.rstrip() for x in sp.communicate()[0].rstrip().split("\n")]
      starts = []
      ends = []
      chrom = sentries[0]['chrom']
      strand = sentries[0]['strand']
      for line in mlines:
        f =  line.rstrip().split("\t")
        tstart = int(f[1])
        tend = int(f[2])
        starts.append(tstart)
        ends.append(tend)
      gline = read + "\t" + read + "\t" + chrom + "\t" + strand + "\t"
      gline += str(starts[0])+"\t" + str(ends[-1])+"\t"
      gline += str(starts[0])+"\t" + str(ends[-1])+"\t"
      gline += ','.join([str(x) for x in starts])+','+"\t"
      gline += ','.join([str(x) for x in ends])+','+"\n"
      glines.append(gline)
    return glines

def getKey(thing):
  return thing['qstart']

def get_locus_group(args,exons):
  #definetly break it up based on chromosome
  chroms = {}
  locus_groups = []
  for entry in exons:
    if entry['chrom'] not in chroms:
      chroms[entry['chrom']] = {}
    if entry['strand'] not in chroms[entry['chrom']]:
      chroms[entry['chrom']][entry['strand']] = []
    chroms[entry['chrom']][entry['strand']].append(entry)
  # now we have entries broken by chromosome
  for chrom in chroms:
   for strand in chroms[chrom]:
    entry_sets = []
    # seed the entries array with all seperate then join on distance
    for entry in chroms[chrom][strand]:
      entry_sets.append([entry])
    new_entries = []
    while len(new_entries) != len(entry_sets):
      #print len(entry_sets)
      if len(new_entries) > 0: entry_sets = new_entries
      new_entries = combine_down(entry_sets)
    entry_sets = new_entries
    #for entry_set in entry_sets:
    #  for entry in entry_set:
    #    print entry
    #  print '----'
    for entry_set in entry_sets:
      locus_groups.append(entry_set)
  return locus_groups

def combine_down(entry_sets):
  ecount = len(entry_sets)
  skip = ecount
  for i in range(0,ecount):
    for j in range(i+1,ecount):
      eset1 = entry_sets[i]
      eset2 = entry_sets[j]
      if near(eset1,eset2,400000):
        new_entries = []
        #print 'near'
        newset = []
        for e in eset1:
          newset.append(e)
        for e in eset2:
          newset.append(e)
        # combine
        for z in range(0,ecount):
          if z == j:
            continue
          elif z == i:
            new_entries.append(newset)
          else:
            new_entries.append(entry_sets[z])
        #  new_entries.append(eset1)
        #print [len(x) for x in new_entries]
        return new_entries
  new_entries = []
  for eset in entry_sets:
    new_entries.append(eset)          
  return new_entries

def near(eset1,eset2,distance):
  for e1 in eset1:
    for e2 in eset2:
      if abs(e1['tstart'] - e2['tstart']) < distance: return True
      if abs(e1['tend'] - e2['tstart']) < distance: return True
      if abs(e1['tstart'] - e2['tend']) < distance: return True
      if abs(e1['tend'] - e2['tend']) < distance: return True
  return False
def make_exons(args,thread_index,thread_count):
  is_sam = True
  if re.search('\.bam$',args.sam_file):
    is_sam = False
  stag = ''
  if is_sam: stag = '-S'
  cmd = 'samtools view -F 4 '+stag+' '+args.sam_file
  spcf = SamBasics.SAMtoPSLconversionFactory()
  if args.reference_genome:
    spcf.set_genome(args.reference_genome)
  sampipe = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
  fname = args.tempdir+'/bedpart.'+str(thread_index)+'.bed'
  of = open(fname,'w')
  z = 0
  with sampipe.stdout as inf:
    for line in inf:
      z+=1
      if z%thread_count != thread_index: continue
      line = line.rstrip()
      if SamBasics.is_header(line):
        continue
      d = SamBasics.sam_line_to_dictionary(line)
      strand = '+'
      if SamBasics.check_flag(d['flag'],16):
        strand = '-'
      seqs = []
      sequence = d['seq']
      seqs.append([d['qname'], d['rname'], strand, d['pos'], d['cigar']])
      m = re.search('XA:Z:(\S+)',line)
      if m and args.use_secondary_alignments:
        e = m.group(1)
        secondaries = e.rstrip(";").split(";")
        for secondary in secondaries:
          m1 = re.match('([^,]+),([+-])(\d+),([^,]+)',secondary)
          if not m1:
            sys.stderr.write("strange secondary format "+secondary+"\n")
            sys.exit()
          seqs.append([d['qname'], m1.group(1),m1.group(2),int(m1.group(3)),m1.group(4)])
      #p.apply_async(get_exons_from_seqs,[seqs,d,spcf])
      exons = get_exons_from_seqs(seqs,d,spcf)
      of.write(exons)
      #return exons
  of.close()

def get_exons_from_seqs(seqs,d,spcf):
  sind = 0
  oline = ''
  for seq in seqs:
    sind+=1
    psec = 'P' #primary or secondary
    if sind > 1: psec = 'S'
    d1 = d.copy()
    d1['rname'] = seq[1]
    if seq[2] == '+':  d1['flag'] = 0
    else: d1['flag'] = 16
    d1['pos'] = seq[3]
    d1['cigar'] = seq[4]
    d1['cigar_array'] = SamBasics.parse_cigar(seq[4])
    skips = set(['H','D','N'])
    total_length = 0
    possible_matches = 0
    indels = 0
    qstart = 0
    if d1['cigar_array'][0]['op'] == 'S':
      qstart = d1['cigar_array'][0]['val']
    if d1['cigar_array'][0]['op'] == 'H':
      qstart = d1['cigar_array'][0]['val']
    for ce in d1['cigar_array']:
      if ce['op'] not in skips:
        total_length += ce['val']
      if ce['op'] == 'M': possible_matches += ce['val']
      elif ce['op'] == 'I':
        indels += ce['val']
      elif ce['op'] == 'D' and ce['val'] < 68:
        indels += ce['val']
    fakeseq = 'N'*total_length
    d1['seq'] = fakeseq
    nline = SamBasics.entry_to_line(d1)
    pline = spcf.convert_line(nline)
    pentry = PSLBasics.line_to_entry(pline)
    #mismatch_count = -1
    #if sind == 1 and args.reference_genome: #for primary alignments we can calculate the number of matches
    #  for i in range(0,len(pentry['blockSizes'])):
    #    tseq = spcf.genome[pentry['tName']][pentry['tStarts'][i]:pentry['tStarts'][i]+pentry['blockSizes'][i]]
    #    qseq = sequence[pentry['qStarts'][i]:pentry['qStarts'][i]+pentry['blockSizes'][i]]
    #    print pentry['blockSizes'][i]
    #    print tseq
    #    print qseq
    #    for j in range(0,len(tseq)):
    #      if tseq[j].upper() != qseq[j].upper(): mismatch_count += 1
    gline = PSLBasics.convert_entry_to_genepred_line(pentry)
    gentry = GenePredBasics.line_to_entry(gline)
    gsmooth = GenePredBasics.smooth_gaps(gentry,68)
    for i in range(0,len(gsmooth['exonStarts'])):
      oline += gsmooth['chrom'] + "\t" + str(gsmooth['exonStarts'][i])+"\t"+str(gsmooth['exonEnds'][i])+"\t"+gsmooth['strand']+"\t"+gsmooth['name']+"\t"+str(possible_matches)+"\t"+str(indels)+"\t"+psec+"\t"+str(qstart)+"\n"
  return oline

if __name__=="__main__":
  main()
