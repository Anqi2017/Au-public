#!/usr/bin/python
import argparse, sys, re, json, multiprocessing
from SequenceBasics import read_fasta_into_hash, GenericFastqFileReader, GenericFastaFileReader, rc
from GenePredBasics import GenePredFile
import PSLBasics

combo_results = []

# This software takes as an input, partial alignments in psl format 
# from programs such as BWA-mem.  If you only have a sam file to begin
# with you can use sam_to_psl.py to get the compatible psl file.
# 
# Hits for a given direction for a given query sequence will be grouped
# together within a given locus range.  For this reason, it is imperative
# that the queries all have unique names.  Therefore the prealignment
# fasta or fastq file will be required and checked for this.

def main():
  parser = argparse.ArgumentParser(description="splice together partial alignments")
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('--fastq_reads')
  group1.add_argument('--fasta_reads')
  parser.add_argument('--genome',help="FASTA reference genome",required=True)
  parser.add_argument('--genepred',help="Transcriptome genepred")
  parser.add_argument('--max_intron_size',type=int,default=100000,help="INT maximum intron size")
  parser.add_argument('--min_intron_size',type=int,default=68,help="INT minimum intron size")
  parser.add_argument('--max_gap_size',type=int,default=10,help="INT gap size in query to join")
  parser.add_argument('--max_search_expand',type=int,default=10,help="INT max search space to expand search for junction")
  parser.add_argument('--direction_specific',action='store_true',help="The direction of the transcript is known and properly oriented already")
  parser.add_argument('--threads',type=int,default=0,help="INT number of threads to use default cpu_count")
  parser.add_argument('-o','--output',default='-',help="FILENAME output results to here rather than STDOUT which is default")
  parser.add_argument('input_alignment',help="FILENAME input .psl file or '-' for STDIN")
  args = parser.parse_args()

  # Read our reference genome
  sys.stderr.write("Reading reference\n")
  ref = read_fasta_into_hash(args.genome)

  # Make sure our reads are unique
  sys.stderr.write("Checking for unqiuely named reads\n")
  reads = check_for_uniquely_named_reads(args) # does a hard exit and error if there are any names repeated
  sys.stderr.write("Reads are uniquely named\n")
  
  # Set number of threads to use
  cpu_count = multiprocessing.cpu_count()
  if args.threads > 0:
    cpu_count = args.threads

  #Set reference splices (if any are available)
  reference_splices = {}
  if args.genepred:
    sys.stderr.write("Reading reference splices from genepred\n")
    reference_splices = get_reference_splices(args)

  sys.stderr.write("Reading alignments into loci\n")

  # Get locus division (first stage)
  # Each read (qName) is separated
  # Then each locus will be specific to at chromosome (tName)
  # Then by (strand), but keep in mind this is the is based on the read
  # Each locus should be specific to a direction but we don't necessarily
  # know direction based on the data we have thus far.  
  inf = sys.stdin
  if args.input_alignment != '-': inf = open(args.input_alignment,'r')
  loci = {}
  for line in inf:
    line = line.rstrip()
    if re.match('^#',line): continue
    psl = PSLBasics.line_to_entry(line)
    if psl['qName'] not in loci:
      loci[psl['qName']] = {}
    if psl['tName'] not in loci[psl['qName']]:
      loci[psl['qName']][psl['tName']] = {}
    if psl['strand'] not in loci[psl['qName']][psl['tName']]:
      loci[psl['qName']][psl['tName']][psl['strand']] = {}
    if psl['tStarts'][0] not in loci[psl['qName']][psl['tName']][psl['strand']]:
      loci[psl['qName']][psl['tName']][psl['strand']][psl['tStarts'][0]] = []
    loci[psl['qName']][psl['tName']][psl['strand']][psl['tStarts'][0]].append(psl)

  sys.stderr.write("breaking loci by genomic distance\n")
  for qname in loci:
    for chr in loci[qname]:
      for strand in loci[qname][chr]:
        #print qname + "\t" + chr + "\t" + strand
        starts = loci[qname][chr][strand].keys()
        current_set = []
        locus_sets = []
        last_end = -1*(args.max_intron_size+2)
        for start in sorted(starts):
          for e in loci[qname][chr][strand][start]:
            start = e['tStarts'][0]+1 # base-1 start of start of alignment
            if start > last_end+args.max_intron_size:
              # we have the start of a new set
              if len(current_set) > 0: 
                locus_sets.append(current_set)
              current_set = []
            last_end = e['tStarts'][len(e['tStarts'])-1]+e['blockSizes'][len(e['tStarts'])-1]
            current_set.append(e)
        if len(current_set) > 0:
          locus_sets.append(current_set)
        loci[qname][chr][strand] = locus_sets # replace what was there with these ordered sets

  locus_total = 0
  for qname in loci:
    for chr in loci[qname]:
      for strand in loci[qname][chr]:
        for locus_set in loci[qname][chr][strand]:
          locus_total+=1  

  sys.stderr.write("Work on each read in each locus with "+str(cpu_count)+" CPUs\n")
  p = multiprocessing.Pool(processes=cpu_count)
  locus_count = 0
  for qname in loci:
    for chr in loci[qname]:
      for strand in loci[qname][chr]:
        #print qname + "\t" + chr + "\t" + strand
        for locus_set in loci[qname][chr][strand]:
          locus_count += 1
          onum = len(locus_set)
          # send blank reference splices unless we have some
          rsplices = {}
          if chr in reference_splices: rsplices = reference_splices[chr]
          #p.apply_async(process_locus_set,(locus_set,args,rsplices,ref[chr],reads[qname],locus_total,locus_count),callback=do_locus_callback)
          r1 = execute_locus(locus_set,args,rsplices,ref[chr],reads[qname],locus_total,locus_count)
          do_locus_callback(r1)
          #nnum = len(new_locus_set)
          #print str(onum) + " to " + str(nnum)
          #for e in new_locus_set:
          #  print PSLBasics.entry_to_line(e)
  p.close()
  p.join() 
  sys.stderr.write("\nfinished\n")

  ofh = sys.stdout
  if not args.output == '-':
    ofh = open(args.output,'w')

  for line in combo_results:
    ofh.write(line)

def execute_locus(locus_set,args,reference_splices,ref,read,locus_total,locus_count):
  # for this locus we will have to check both directions
  print "before: " + str(len(locus_set))
  print [int(x['matches']) for x in locus_set]
  print max([int(x['matches']) for x in locus_set])
  rsplices_plus = {}
  rsplices_minus = {}
  if '+' in reference_splices: rsplices_plus = reference_splices['+']
  if '-' in reference_splices: rsplices_minus = reference_splices['-']
  [r1 ,lt1, lc1] = process_locus_set(locus_set,args,rsplices_plus,ref,read,locus_total,locus_count,'+')
  [r2, rt2, rc2] = process_locus_set(locus_set,args,rsplices_minus, ref, read, locus_total, locus_count, '-')
  print "after: " 
  print "  positive"
  print "  "+str([int(x['matches']) for x in r1])
  print "  "+str([int(x['misMatches']) for x in r1])
  print "  "+str([int(x['tBaseInsert']) for x in r1])
  print "  "+str([int(x['qBaseInsert']) for x in r1])
  print "  negative"
  print "  "+str([int(x['matches']) for x in r2])
  print "  "+str([int(x['misMatches']) for x in r2])
  print "  "+str([int(x['tBaseInsert']) for x in r2])
  print "  "+str([int(x['qBaseInsert']) for x in r2])

  #print r1
  #print r2
def do_locus_callback(cbr):
  if not cbr: return
  [r,tot,cnt] = cbr
  global combo_results
  for e in r:
    combo_results.append(PSLBasics.entry_to_line(e)+"\n")
  sys.stderr.write(str(cnt)+"/"+str(tot)+"  \r")
  return

# now that we have our best choice, lets put them together
def do_combine_operation(best_option,left,right,read,seq,args):
  #print "choice is "+str(best_option)
  left_target = best_option[0]
  right_target = best_option[1]
  left_query = best_option[2]
  right_query = best_option[3]
  # store for output
  q_start_array = []
  t_start_array = []
  block_size_array = []

  left_query_start = left['qStarts'][0]
  left_target_start = left['tStarts'][0]
  for i in range(0,len(left['tStarts'])):
    tstart = left['tStarts'][i]
    tend = left['tStarts'][i]+left['blockSizes'][i]
    qstart = left['qStarts'][i]
    qend = left['qStarts'][i]+left['blockSizes'][i]
    if left_query <= qstart+1: break
    left_query_start = qstart
    left_target_start = tstart
    if left_query <= qend: break
    q_start_array.append(qstart)
    t_start_array.append(tstart)
    block_size_array.append(left['blockSizes'][i])

  #print "left things"
  #print [left_query_start+1,left_query]
  #print [left_target_start+1,left_target]

  right_query_end = right['qStarts'][0]+right['blockSizes'][0]
  right_target_end = right['tStarts'][0]+right['blockSizes'][0]
  right_outer_index = 0
  for j in range(0,len(right['tStarts'])):
    tstart = right['tStarts'][j]
    tend = right['tStarts'][j]+right['blockSizes'][j]
    qstart = right['qStarts'][j]
    qend = right['qStarts'][j]+right['blockSizes'][j]
    right_outer_index = j+1
    if right_query <= qstart+1: break
    right_query_end = qend
    right_target_end = tend
    if right_query < qend: break
  #print "right things"
  #print [right_query+1,right_query_end]
  #print [right_target+1,right_target_end]
  working_read = read.upper()
  if left['strand'] == '-': working_read = rc(read.upper())
  pread = working_read[left_query_start:right_query_end]
  tseq = seq[left_target_start:left_target].upper()+seq[right_target-1:right_target_end].upper()
  res = needleman_wunsch(pread,tseq)
  #print "short needleman wunsch"
  #print res[0]
  #print res[1]

  # Fun part of making the new portion of the alignment
  qindex = left_query_start
  tindex = left_target_start
  in_alignment = 0
  alignment = None
  bynumbers = None
  for i in range(0,len(res[0])):
    if res[0][i] == '-':  #insertion in target (gap in query)
      tindex += 1
      in_alignment = 0
    elif res[1][i] == '-':  #insertion in query (gap in target)
      qindex += 1
      in_alignment = 0
    else: # we are in an alignment
      if in_alignment == 0:
        # output buffered result
        if alignment:
          if len(alignment[0]) > 0:
            q_start_array.append(bynumbers[0])
            t_start_array.append(bynumbers[1])
            block_size_array.append(len(alignment[0]))
        alignment = ['','']
        bynumbers = [qindex,tindex,qindex,tindex]
      in_alignment = 1
      alignment[0] += res[0][i]
      alignment[1] += res[1][i]
      bynumbers[2] += 1
      bynumbers[3] += 1
      qindex+=1
      tindex+=1
    if qindex == right_query: # switch forward
      #print "switch"
      #print str(tindex) + "\t" + str(right_target)
      #print str(qindex) + "\t" + str(right_query)
      if not tindex == right_target: 
        in_alignment = 0
      tindex = right_target
  if alignment:
    if len(alignment[0]) > 0:
      q_start_array.append(bynumbers[0])
      t_start_array.append(bynumbers[1])
      block_size_array.append(len(alignment[0]))
    #print bynumbers


  for i in range(right_outer_index,len(right['blockSizes'])):
    q_start_array.append(right['qStarts'][i])
    t_start_array.append(right['tStarts'][i])
    block_size_array.append(right['blockSizes'][i])

  #now we can finally construct a psl line
  #we won't keep track of repeats for now
  matches = 0
  misMatches = 0
  repMatches = 0
  nCount = 0
  qNumInsert = 0
  qBaseInsert = 0
  tNumInsert = 0
  tBaseInsert = 0
  strand = left['strand']
  qName = left['qName']
  qSize = len(read)
  qStart = q_start_array[0]
  qEnd = q_start_array[len(q_start_array)-1]+block_size_array[len(block_size_array)-1]
  tName = left['tName']
  tSize = len(seq)
  tStart = t_start_array[0]
  tEnd = t_start_array[len(t_start_array)-1]+block_size_array[len(block_size_array)-1]
  blockCount = len(block_size_array)
  blockSizes = ','.join([str(x) for x in block_size_array])+','
  qStarts = ','.join([str(x) for x in q_start_array])+','
  tStarts = ','.join([str(x) for x in t_start_array])+','

  prev_q_end = None
  prev_t_end = None
  for i in range(0,len(block_size_array)):
    qseg = working_read[q_start_array[i]:q_start_array[i]+block_size_array[i]]
    tseg = seq[t_start_array[i]:t_start_array[i]+block_size_array[i]].upper()
    for j in range(0,len(qseg)):
      if qseg[j] == 'N': nCount += 1
      if qseg[j] == tseg[j]: matches += 1
      else:
        misMatches += 1
    if prev_t_end:
      t_dist = t_start_array[i]-prev_t_end
      if t_dist > 0 and t_dist < args.min_intron_size: #we have an insert into the target and its not an intron
        tNumInsert += 1
        tBaseInsert += t_dist
    if prev_q_end:
      q_dist = q_start_array[i]-prev_q_end
      if q_dist > 0:
        qNumInsert += 1
        qBaseInsert += q_dist
    prev_q_end = q_start_array[i]+block_size_array[i]
    prev_t_end = t_start_array[i]+block_size_array[i]

  # now we have everything to make the line
  combo_line = str(matches) + "\t" + str(misMatches) + "\t" + str(repMatches) + "\t" \
             + str(nCount) + "\t" + str(qNumInsert) + "\t" + str(qBaseInsert) + "\t" \
             + str(tNumInsert) + "\t" + str(tBaseInsert) + "\t" \
             + strand + "\t" + qName + "\t" + str(qSize) + "\t" \
             + str(qStart) + "\t" + str(qEnd) + "\t" \
             + tName + "\t" + str(tSize) + "\t" \
             + str(tStart) + "\t" + str(tEnd) + "\t" + str(blockCount) + "\t" \
             + blockSizes + "\t" + qStarts + "\t" + tStarts
  #print combo_line
  #print q_start_array
  #print t_start_array
  #print block_size_array
  #  print str(right['qStarts'][i])+"\t"+str(right['qStarts'][i]+right['blockSizes'][i])
  #  print i
  return PSLBasics.line_to_entry(combo_line)

# Here is the heart of the program
def process_locus_set(locus_set,args,reference_splices,seq,read,locus_total,locus_count,orientation):
  print orientation
  lcount = len(locus_set)
  #print len(reference_splices)
  #print '----'
  #print lcount
  score_set = []
  if lcount == 1: # only one entry so nothing to combine
    return [locus_set, locus_total, locus_count]
  #for speed lets do greedy joining of alignments
  stayin = True
  while stayin == True:
    stayin = False
    lcount = len(locus_set)
    buffer = []
    for i in range(1,lcount):
      left = locus_set[i-1].copy()
      right = locus_set[i].copy()
      #if orientation == '+':   # do it like normal
      #  combo = combine(left,right,args,reference_splices,seq,read,orientation)
      #else:  # then its the other orientation
      #  nstrand = opposite(left['strand'])
      #  left['strand'] = nstrand
      #  right['strand'] = nstrand
      #  nread = rc(read)
      #  combo = combine(left,right,args,reference_splices,seq,read,orientation)
      combo = combine(left,right,args,reference_splices,seq,read,orientation)
      if not combo:
        buffer.append(left)
        #print "couldn't combine"
        continue
      else:
        #print "could combine" 
        stayin = True  # continue looping if a change is made
      newset = buffer[:]
      newset.append(combo)
      for j in range(i+1,lcount):
        newset.append(locus_set[j])
      locus_set = newset
      break
  return [locus_set, locus_total, locus_count]
  #print len(locus_set)
  #print '---'

def opposite(strand):
  if strand == '+': return '-'
  return '+'

def combine(left,right,args,reference_splices,seq,read,orientation):
  #Kind of the business end where we combine two psl entries
  # Perform a check for an overlap (or near gap) sufficient for consideration
  if left['qEnd'] < right['qStart']-args.max_gap_size: # no overlap may want to have a seperate parameter for a max gap size
    return None
  target_options = get_options(left, \
                   min(left['qEnd'],right['qStart']+1), \
                   left['qEnd'], \
                   right, \
                   right['qStart']+1, \
                   max(left['qEnd'],right['qStart']+1))                     


  # Try to find the junction site
  junction_choice = None
  reference_options = []
  for j in [json.loads(x) for x in reference_splices]:
    if j[0] in target_options:
      if j[1] in target_options[j[0]]:
        for op in target_options[j[0]][j[1]]:
          dist = abs(op[2]-1)+op[3]+op[4]
          # this is where we keep ourselves from looking too far away
          if dist < args.max_search_expand:
            reference_options.append([j[0],j[1]] + op)
  # check for a cannonical spice site
  strand = left['strand']
  candidate_options = []
  for l_t in target_options:
    for r_t in target_options[l_t]:
      candidate = seq[l_t:l_t+2].upper()+'-'+seq[r_t-3:r_t-1].upper()
      iscan = False
      if strand == '+' and is_canon(candidate):
        iscan = True
      if strand == '-' and is_revcanon(candidate):
        iscan = True
      if iscan:
        for entry in target_options[l_t][r_t]:
          dist = abs(entry[2]-1)+entry[3]+entry[4]
          # this is where we keep ourselves from looking too far away
          if dist < args.max_search_expand:
            candidate_options.append([l_t,r_t]+entry+[strand,candidate])

  #for c in candidate_options: print c
  if len(candidate_options) == 0: return None

  # For choosing the best candidate we don't need to align all of the seq
  prefered_alignment_length = 50
  # Which left alignment segment has candidates in it?
  leftlen = len(left['qStarts'])
  minleft_query = min([x[2] for x in candidate_options])
  left_nearest = 0
  for i in range(0,leftlen):
    if minleft_query >= left['qStarts'][i] + 1:
      left_nearest = i
    else:
      break
  
  left_choice = left_nearest
  if left_nearest > 0:
    for i in range(left_nearest,0-1,-1):
      tot = left['qStarts'][left_nearest]-left['qStarts'][i]
      if tot > prefered_alignment_length:
        left_choice = i
        break
      left_choice = i

  # Which right alignment segment has candidates in it?
  rightlen = len(right['qStarts'])
  maxright_query = max([x[3] for x in candidate_options])
  right_nearest = 0
  for i in range(0,rightlen):
    if maxright_query <= right['qStarts'][i]+right['blockSizes'][i]:
      break
    else:
      right_nearest = i


  right_choice = right_nearest
  if right_nearest > 0:
    for i in range(right_nearest,rightlen):
      tot = (right['qStarts'][i]+right['blockSizes'][i])-(right['qStarts'][right_nearest]+right['blockSizes'][right_nearest])
      if tot > prefered_alignment_length:
        right_choice = i
        break
      right_choice = i

  # now left_choice and right_choice contain bounds to align

  # we can come up with options based on a needleman wunsch across the entire thing
  wread = read
  if strand == '-': wread = rc(read)
  newoptions = []
  for option in reference_options + candidate_options:
    #leftlen = len(left['qStarts'])
    #rightlen = len(right['qStarts'])
    [a1,a2,score] = needleman_wunsch( \
                    wread[left['qStarts'][left_choice]:right['qStarts'][right_choice]+right['blockSizes'][right_choice]].upper(), \
                    seq[left['tStarts'][left_choice]:option[0]].upper() + \
                    seq[option[1]-1:right['tStarts'][right_choice]+right['blockSizes'][right_choice]].upper())
    newoptions.append([score,a1,a2]+[option])

  #best_option = get_best_option(reference_options, candidate_options, args, strand)
  #if not best_option:
  #  return None
  [best_option,best_score,best_align] = get_best_option2(newoptions)
  if not best_option:
    return None
  #print best_option
  #print best_align[0]
  #print best_align[1]
  combo = do_combine_operation(best_option,left,right,read,seq,args)
  return combo

def get_best_option2(newoptions):
  bestscore = 0
  bestop = None
  besta = None
  for [score,a1,a2,op] in newoptions:
    if score > bestscore:
      bestscore = score
      bestop = op
      besta = [a1,a2]
  return bestop,bestscore,besta

# This function could probably use a project of its own
# but we will be greedy for now and just take the first candidate
# with the smallest gap in the query alignments, that agrees with
# first either a reference junction
# or second a canonincal splice signal (starting with the major signal)
# if the major signal is not detected in the search space, we expand to
# the minor signals.  
def get_best_option(refop,canop,args,strand):
  canon_set = ['GT-AG','GC-AG','AT-AC'] # order by priority
  noncanon_set = ['CT-AC','CT-GC','GT-AT']
  use_set = canon_set
  if strand == '-': use_set = noncanon_set
  for canon in use_set:
    for i in range(0,args.max_gap_size+1):
      for op in refop:
        distance = abs(1-op[4])+op[5]+op[6] #distance from optimum
        #print distance
        #print op
        if distance == i:
          return op
      for op in canop:
        if not op[8] == canon: continue # only do the best site at a time
        distance = abs(1-op[4])+op[5]+op[6]
        #print distance
        #print op
        if distance == i:
           return op
  return None

def is_canon(input):
  v = set()
  v.add('GT-AG')
  v.add('GC-AG')
  v.add('AT-AC')
  if input in v: return True
  return False

def is_revcanon(input):
  v = set()
  v.add('CT-AC')
  v.add('CT-GC')
  v.add('GT-AT')
  if input in v: return True
  return False


def is_cannonical(chr,strand,left,right,seq):
  
  return True

# Get mapping between target and query, for aligned bases only
def get_basic_options(psl,qmin,qmax):
  options = {}
  query_options = {}
  for i in range(0,len(psl['blockSizes'])):
    for j in range(0,psl['blockSizes'][i]):
      q  = psl['qStarts'][i]+j+1
      t  = psl['tStarts'][i]+j+1
      if q >= qmin and q <= qmax:  
        options[t] = q
        query_options[q] = t
  return [options, query_options]


#Between qmin and qmax (1-base query index) for a psl entry
#List all target coordinates
def get_options(left,lqmin,lqmax,right,rqmin,rqmax):
  left_options = {}
  left_query_options = {}

  inner_expand = 10
  outer_expand = 10
  # get basic mapping
  [left_options, left_query_options] = get_basic_options(left,lqmin-inner_expand,lqmax)
  right_options = {}
  right_query_options = {}
  [right_options, right_query_options] = get_basic_options(right,rqmin,rqmax+inner_expand)

  # now we want to expand this somewhat.  
  # For all the points of interest in the query, we want to
  # know the distance to the overlapping region
  # know the target sequences associated with that point 
  # This part will force an assignment to a query and target for all bases
  # in our search space
  [left_full_target_options, left_full_query_options, ltmin, ltmax] = expand_options(left_options, left_query_options,lqmin-inner_expand,lqmax)
  [right_full_target_options, right_full_query_options, rtmin, rtmax] = expand_options(right_options, right_query_options,rqmin,rqmax+inner_expand)

  [left_expanded_target_options, left_expanded_query_options] = expand_gap(left_full_target_options, left_full_query_options,ltmax, ltmax+1, ltmax+outer_expand)
  [right_expanded_target_options, right_expanded_query_options] = expand_gap(right_full_target_options, right_full_query_options,rtmin, rtmin-outer_expand-1, rtmin-1)

  assessment = {}
  for l_t in left_expanded_target_options:
    for r_t in right_expanded_target_options:
      for l_q in left_expanded_target_options[l_t]:
        for r_q in right_expanded_target_options[r_t]:
          distance = r_q-l_q
          if l_t not in assessment:
            assessment[l_t] = {}
          if r_t not in assessment[l_t]:
            assessment[l_t][r_t] = []
          lgap = 0
          if l_t > ltmax: lgap = l_t-ltmax
          rgap = 0
          if r_t < rtmin: rgap = rtmin-r_t
          assessment[l_t][r_t].append([l_q,r_q,distance,lgap,rgap])
  return assessment

def expand_options(options, query_options,qmin,qmax):
  previous_target = None
  full_query_options = {}
  for q in range(qmin,qmax+1):
    target = previous_target
    if q in query_options:
      target = query_options[q]
      previous_target = target
    if q not in full_query_options:
      full_query_options[q] = []
    if target: full_query_options[q].append(target)

  full_target_options = {}
  tmin = min(options.keys())
  tmax = max(options.keys())
  previous_query = None
  for t in range(tmin,tmax+1):
    query = previous_query
    if t in options:
      query = options[t]
      previous_query = query
    if t not in full_target_options:
      full_target_options[t] = []
    if query: full_target_options[t].append(query)
  return [full_target_options, full_query_options, tmin, tmax]


# now we can expand by getting the max target
def expand_gap(full_target_options, full_query_options,border,r1,r2):
  queries = full_target_options[border]
  for i in range(r1,r2+1):
    for q in queries: #probably only one
      full_query_options[q].append(i)
      if i not in full_target_options:
        full_target_options[i] = []
      full_target_options[i].append(q)
  return [full_target_options, full_query_options]

def get_reference_splices(args):
    reference_splices = {}
    gpf = GenePredFile(args.genepred)
    for ge in gpf.entries:
      chrom = ge.entry['chrom']
      strand = ge.entry['strand']
      junctions = ge.junctions
      for junc in junctions:
        m = re.match('^(\S+):(\d+),(\S+):(\d+)$',junc)
        if not m: continue
        if chrom not in reference_splices:
          reference_splices[chrom] = {}
        if strand not in reference_splices[chrom]:
          reference_splices[chrom][strand] = set()
        left_coord = int(m.group(2))
        right_coord = int(m.group(4))
        reference_splices[chrom][strand].add(json.dumps([left_coord,right_coord]))
    return reference_splices

# exits execution of the program if there are repeat names of reads
# returns a hash of reads
def check_for_uniquely_named_reads(args):
  observed_reads = set()
  reads = {}
  if args.fastq_reads:
    gfr = GenericFastqFileReader(args.fastq_reads)
    while True:
      e = gfr.read_entry()
      if not e: break
      reads[e['name']] = e['seq']
      if e['name'] in observed_reads:
        sys.stderr.write("ERROR observed reads must be uniquely named")
        sys.exit()
      observed_reads.add(e['name'])
  elif args.fasta_reads:
    gfr = GenericFastaFileReader(args.fasta_reads)
    while True:
      e = gfr.read_entry()
      if not e: break
      reads[e['name']] = e['seq']
      if e['name'] in observed_reads:
        sys.stderr.write("ERROR observed reads must be uniquely named")
        sys.exit()
      observed_reads.add(e['name'])
  return reads


def needleman_wunsch(s1,s2):
  F = []
  d = -15
  for i in range(0,len(s1)+1):
    temp = []
    F.append(temp)
    for j in range(0,len(s2)+1):
      F[i].append(0)
  for i in range(0,len(s1)+1):
    F[i][0] = d*i
  for j in range(0,len(s2)+1):
    F[0][j] = d*j
  for i in range(1,len(F)):
    for j in range(1,len(F[i])):
      match = F[i-1][j-1]+Snw(s1[i-1],s2[j-1])
      deletion = F[i-1][j]+d
      insertion = F[i][j-1]+d
      F[i][j] = max(match,deletion,insertion)
  a1 = ''
  a2 = ''
  i = len(s1)
  j = len(s2)
  while i > 0 or j > 0:
    if i > 0 and j > 0 and F[i][j] == F[i-1][j-1]+Snw(s1[i-1],s2[j-1]):
      a1 = s1[i-1] + a1
      a2 = s2[j-1] + a2
      i -= 1
      j -= 1
    elif i > 0 and F[i][j] == F[i-1][j] + d:
      a1 = s1[i-1] + a1
      a2 = '-' + a2
      i -= 1
    else:
      a1 = "-" + a1
      a2 = s2[j-1]+a2
      j -= 1
  return [a1,a2,F[len(s1)][len(s2)]] 

def Snw(c1,c2):
  if c1 == c2: return 10
  else: return -5

if __name__=="__main__":
  main()
