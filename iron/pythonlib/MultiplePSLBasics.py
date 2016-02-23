import re, sys
from RangeBasics import GenomicRange, Bed
from SequenceBasics import rc as rc_seq
from PSLBasics import is_valid, PSL
import GraphBasics

class MultiplePSLAlignments:
  def __init__(self):
    self.entries = [] # a list of PSL entries
    self.qName = None
    self.g = None
  def entry_count(self):
    return len(self.entries)
  def add_entry(self,entry):
    if not self.qName:
      self.qName = entry.value('qName')
    else:
      if entry.value('qName') != self.qName:
        sys.stderr.write("WARNING multiple alignments must have the same query name.  This entry will not be added\n")
        return False
    self.entries.append(entry)
    return True
  def get_root_paths(self):
    if not self.g: 
      sys.stderr.write("ERROR must run compatible_graph() before get_root_paths()\n")
      sys.exit()
    rs = self.g.get_roots()
    paths = []
    for r in rs:
      ps = self.g.get_directed_paths_from_node(r)
      for p in ps:
        paths.append(p)
    outputs = []
    for path in paths:
      simples = _traverse([self.g.get_nodes()[x].get_payload() for x in path])
      for simple in simples:
        #print str(simple)+" "+str([self.entries[i].get_target_bed().get_range_string() for i in simple])+" "+str([self.entries[i].get_query_bed().get_range_string() for i in simple])
        outputs.append(simple)
    return outputs

  # Pre: max_intron - set this to -1 if you don't care about target placement (use for fusion)
  #      max_gap - set to -1 to allow unaligned sequences of any size on the query
  #      max_query_overlap - The largest number of bases you permit two query sequences to overlap
  #      max_target_overlap - The largest number of bases you permit two target sequences to overlap
  def compatible_graph(self,max_intron=400000,max_gap=-1,max_query_overlap=0,max_target_overlap=0,max_query_fraction_overlap=-1):
    # Create a flow graph of how multiple PSL files can be connected
    g = GraphBasics.Graph() #directed graph
    node_dict = {}
    for i in range(0,len(self.entries)):
      n1 = GraphBasics.Node([i])
      g.add_node(n1)
      node_dict[i] = n1
    for i in range(0,len(self.entries)):
      for j in range(0,len(self.entries)):
        if i == j: continue # not looking when its itself
        #See if we can go from i to j
        qd = self.entries[i].query_distance(self.entries[j])
        td = self.entries[i].target_distance(self.entries[j],use_direction=True)
        if max_intron >= 0 and td == -1: continue
        if max_intron > 0:
          if td > max_intron: continue
        if max_gap > 0:
          if qd > max_gap: continue
        target_overlap = self.entries[i].target_overlap_size(self.entries[j])
        query_overlap = self.entries[i].query_overlap_size(self.entries[j])
        if max_query_fraction_overlap > 0:
          frac = max(float(query_overlap)/float(self.entries[i].get_coverage()),float(query_overlap)/float(self.entries[j].get_coverage()))
          if frac > max_query_fraction_overlap: continue
        if query_overlap > max_query_overlap and max_query_overlap >= 0: continue
        if target_overlap > max_target_overlap and max_target_overlap >= 0: continue
        # make sure j target is greater than i target
        # check for plus stand
        if self.entries[j].value('tStarts')[0] < self.entries[i].value('tStarts')[0] and max_intron >= 0 and self.entries[i].value('strand')=='+':
          continue
        # check for minus strand
        if self.entries[j].value('tStarts')[0] > self.entries[i].value('tStarts')[0] and max_intron >= 0 and self.entries[i].value('strand')=='-':
          continue
        # make sure j query is greater than i query
        if self.entries[j].value('qStarts_actual')[0] < self.entries[i].value('qStarts_actual')[0]:
          continue
        #print self.entries[i].value('tName') +"\t"+self.entries[j].value('tName')+"\t"+str(qd)+"\t"+str(td)
        # add to graph
        #n1 = GraphBasics.Node([i])
        #n2 = GraphBasics.Node([j])
        g.add_edge(GraphBasics.Edge(node_dict[i],node_dict[j]))
    g.merge_cycles()
    self.g = g
    return g

# values is an array of possible values at each positoin
# return possible arrays
def _traverse(values,start_index=0,prev=[]):
  if start_index >= len(values): 
    return [prev]
  outputs = []
  for v in values[start_index]:
    newprev = prev[:]
    newprev.append(v)
    output = _traverse(values,start_index=start_index+1,prev=newprev)
    if output:
     for o in output:
      outputs.append(o)
  return outputs

class MultiplePSLAlignmentsOld:
  def __init__(self):
    self.entries = [] # a list of PSL entries
    self.minimum_coverage = 1 #how many base pairs an alignment must cover to be part of a multiple alignment
    self.qName = None
    self.best_coverage_fraction = 0.9 #how much of an alignment be the best alignment
                                       #where it aligns to be kept later on
    self.multiply_mapped_minimum_quality_difference = 0.05 # if the qualities are similar we can
                                                           # treat it as an ambiguous alignment
    self.multiply_mapped_minimum_overlap = 50 # The minimum number of bases that must be shared between two alignments to consider them multiply mapped
    self.best_alignment = None # This will hold a BestAlignmentCollection class and is set by 'best_query'
    self.verbose = False
  def entry_count(self):
    return len(self.entries)
  def set_minimum_coverage(self,mincov):
    self.minimum_coverage = mincov
  def add_line(self,line):
    self.add_entry(line_to_entry(line))
  def add_entry(self,entry):
    if not self.qName:
      self.qName = entry.value('qName')
    else:
      if entry.value('qName') != self.qName:
        sys.stderr.write("WARNING multiple alignments must have the same query name.  This entry will not be added\n")
        return False
    if self.minimum_coverage > 1:
      cov = entry.get_coverage()
      if cov < self.minimum_coverage:
        if self.verbose: sys.stderr.write("WARNING alignment less than minimum coverage.\n")
        return False
    self.entries.append(entry)
    return True
  def get_alignment_count(self):
    return len(self.entries)
  def get_tNames(self):
    names = set()
    for name in [x.value('tName') for x in self.entries]:
      names.add(name)
    return sorted(list(names))

  # Use the multiple alignments to set information about the query
  # Pre: alignment(s) to have been loaded alread
  # Post: A hash by position of query that contains all
  #       alignments and information on the quality of the alignment
  #       self.query[query_base index-0][entry index]
  #       then contains tName, coverage, and quality
  def populate_query(self):
    query = {}
    # Go through each alignment
    for eindex in range(0,len(self.entries)): 
      e = self.entries[eindex].entry
      cov = self.entries[eindex].get_coverage()
      qual = self.entries[eindex].get_quality()
      # Go through each block of the alignment
      for i in range(0,e['blockCount']):
        # Set relevant mapped alignments for each base of the query
        for z in range(e['qStarts_actual'][i],e['qStarts_actual'][i]+e['blockSizes'][i]):
          if z not in query: query[z] = {}
          query[z][eindex] = {}
          query[z][eindex]['tName'] = e['tName']
          query[z][eindex]['coverage'] = cov
          query[z][eindex]['quality'] = qual
    self.query = query
    return

  # Try to determine the fraction of each alignment that is 'best' for
  # Wherever it is aligning to.  
  # Pre: 1. A list of indecies of valid alignment entries in self.entries
  #      2. A list of indecies valid bases in self.quality 
  # Post: A hash keyed by valid entry indecies that contains
  #       'coverage', 'quality' and 'best_coverage'
  def evaluate_entry_coverage(self,valid_entries,valid_bases):
    qbases = sorted(valid_bases)
    es = {}
    for i in valid_entries:
      es[i]={}
      es[i]['coverage'] = self.entries[i].get_coverage()      
      es[i]['quality'] = self.entries[i].get_quality()      
      es[i]['best_coverage'] = 0 # need to calculate this
    bestbases = {}
    # Step 1:  Calculate coverage fraction for all alignments
    for z in qbases:
      bestindex = -1
      bestquality = 0
      for eindex in self.query[z]:
        if eindex not in valid_entries: continue # only consider the candidates
        if es[eindex]['quality'] > bestquality:
          bestindex = eindex
          bestquality = self.query[z][eindex]['quality']
      if bestindex > -1:
        bestbases[z] = bestindex
    # For each alignment set the amount that alignment constitutes the best 
    # alignment
    for z in sorted(bestbases.keys()):
      ebest = bestbases[z]
      es[ebest]['best_coverage'] += 1
    return [es,bestbases]    

  # Filter our best based on a requirement for being the 'best' for some large
  # fraction of an aligned region.  Reevaluates best coverage for remaining
  # Pre: 1. list of valid entries
  #      2. list of valid bases
  # Post: entry evaluation for entries after removing filtered entries
  #       bases post filtering
  def filter_by_coverage_fraction(self,valid_entries,qbases):
    filteredbases = {}
    contributing_indecies = set()
    [entry_evaluation,temp] = self.evaluate_entry_coverage(valid_entries,qbases)
    # Step 2: Filter out alignments not satisfying the coverage fraction
    for z in sorted(qbases):
      bestindex = -1
      bestquality = 0
      for eindex in entry_evaluation.keys():
        if eindex not in self.query[z]: continue
        if float(entry_evaluation[eindex]['best_coverage'])/float(entry_evaluation[eindex]['coverage']) < self.best_coverage_fraction: continue
        if entry_evaluation[eindex]['quality'] > bestquality:
          bestindex = eindex
          bestquality = self.query[z][eindex]['quality']
      if bestindex > -1:
        filteredbases[z] = bestindex
        contributing_indecies.add(bestindex)
    nentries = list(contributing_indecies)
    [new_eval,new_bases] = self.evaluate_entry_coverage(nentries,filteredbases.keys())
    return [new_eval,new_bases]
  # Very simply return one psl entry with the best coverage
  # Pre: entry has been added
  # Post: a PSL type
  def best_psl(self):
    best = 0
    for e in self.entries:
      if e.value('matches') > best: best = e.value('matches')
    for e in self.entries:
      if e.value('matches') == best:
        return e.copy()

  # Read throught he query data and find the best explainations
  # Pre: 1.  Have loaded in alignments
  #      2.  Have ran populate_query()
  #      Now populate query contains a hash of query indecies
  #      that have the hashes of matching alignments
  #      and information regarding the quality of each of those alignments
  # Post: Sets self.best_contributing_entries 
  #        and self.best_alignment_segments = []
  def best_query(self):
    qbases = sorted(self.query.keys())
    all_entries = range(0,len(self.entries))
    # Step 1:  Calculate coverage fraction for all alignments
    [entry_evaluation,bases] = self.evaluate_entry_coverage(all_entries,qbases)
    # Step 2: Filter out alignments not satisfying the coverage fraction
    [filtered_entry_evaluation,filtered_bases] = self.filter_by_coverage_fraction(entry_evaluation,bases.keys())

    # Get bed information for the alignment
    qbases = sorted(filtered_bases.keys())
    if len(qbases) == 0: return False
    qstart = qbases[0]
    current = filtered_bases[qstart]
    last = qstart
    beds = []
    eindex = filtered_bases[qstart]
    for i in range(1,len(qbases)):
      e = self.entries[eindex]
      current = filtered_bases[qbases[i]]
      if current not in filtered_entry_evaluation: continue
      if eindex != current:
        beds.append([qstart,last+1,eindex])
        qstart = qbases[i]
      eindex = filtered_bases[qbases[i]]
      last = qbases[i]
    beds.append([qstart,last+1,eindex])
    contributing_indecies = set()
    filtered_beds = []
    for bed in beds:
      seglen = bed[1]-bed[0]
      if seglen < self.minimum_coverage: # remove a fragment too short to call
        for z in range(bed[0],bed[1]):
          if z in filtered_bases:
            del filtered_bases[z]
      else:
        filtered_beds.append(bed)
        contributing_indecies.add(bed[2])
    #if len(contributing_indecies) < 2: return
    #print '---'
    #for i in sorted(list(contributing_indecies)):
    #  print str(i)+"\t"+self.entries[i]['tName']+"\t"+self.entries[i]['strand']+"\t"+str(get_coverage(self.entries[i]))+"\t"+str(get_quality(self.entries[i]))
    self.best_alignment = BestAlignmentCollection()
    entries_present = set()
    for bed in filtered_beds:
      temp = {}
      temp['query_bed'] = [bed[0],bed[1]]
      temp['psl_entry_index'] = bed[2]
      temp['multiply_mapped'] = self.get_multiply_mapped(bed[2],bed[0],bed[1])
      entries_present.add(bed[2])
      self.best_alignment.segments.append(temp)
    for i in entries_present:
      self.best_alignment.entries[i] = self.entries[i]
    self.best_alignment.qName = self.qName
    return self.best_alignment

  def get_multiply_mapped(self,eindex,bed_start,bed_finish):
    multibase = {}
    for i in range(bed_start,bed_finish):
      if i in self.query:
        if eindex in self.query[i]:
          bestquality = self.query[i][eindex]['quality']
          for eindex2 in self.query[i]:
            if eindex2 == eindex: continue
            if eindex2 not in multibase: multibase[eindex2] = 0
            if self.query[i][eindex]['quality'] > bestquality - self.multiply_mapped_minimum_quality_difference:
              multibase[eindex2] += 1
            if multibase[eindex2] >= self.multiply_mapped_minimum_overlap:
              return True
    return False

def get_psl_quality(entry):
  return float(entry['matches'])/float(entry['matches']+entry['misMatches']+entry['tNumInsert']+entry['qNumInsert'])

# Store the result of a 'best_query' in this
# Can go on to calculate get_trimmed_entries() to cut our entries down by segment
class BestAlignmentCollection:
  def __init__(self):
    self.entries = {}  # psl entries stored by an integer key
    self.segments = [] # contains a query_bed and a psl_entry_index
    self.qName = None
    self.minimum_overlap = 1 # by default consider any overlap as reportable overlap
    self.overlapping_segment_targets = None # set by find_overlapping_segment_targets
    self.minimum_locus_distance = 400000 # minimum number of bases to consider something a different locus 
    self.segment_trimmed_entires = None # set by function can be set to an array equal to size segments
    return

  # Pre: A best alignment collection, for each segment, trim the PSL entry
  #      to fit within these query bed bounds
  # Post: sets self.segement_trimmed_entries
  def get_trimmed_entries(self):
    self.segment_trimmed_entries = []
    for seg in self.segments:
      qbed = seg['query_bed']
      psl = self.entries[seg['psl_entry_index']]
      tpsl = psl.left_qactual_trim(qbed[0]+1)
      tpsl = tpsl.right_qactual_trim(qbed[1])
      self.segment_trimmed_entries.append(tpsl)
    return self.segment_trimmed_entries

  def has_multiply_mapped_segments(self):
    for i in range(0,len(self.segments)):
      if self.segments[i]['multiply_mapped']: return True
    return False

  def has_overlapped_segments(self):
    if not self.overlapping_segment_targets:
      self.find_overlapping_segment_targets()
    if len(self.overlapping_segment_targets.keys()) > 0:
      return True
    return False

  def segment_count(self):
    return len(self.segments)

  def alignment_count(self):
    return len(self.entries)

  def get_gap_sizes(self):
    if len(self.segments)==0: return [0]
    return [self.segments[x]['query_bed'][0]-self.segments[x-1]['query_bed'][1] for x in range(1,len(self.segments))]
  def print_report(self):
    if not self.overlapping_segment_targets:
      self.find_overlapping_segment_targets()
    print '-----'
    print self.qName
    if len(self.entries) > 1:
      biggest_gap_between_entries = max(self.get_gap_sizes())
      print str(biggest_gap_between_entries)+" biggest gap between entries"
    for i in range(0,len(self.segments)):
      overstring = ''
      if i in self.overlapping_segment_targets: overstring = 'OVERLAPPED'
      eindex = self.segments[i]['psl_entry_index']
      mm = self.segments[i]['multiply_mapped']
      mmstring = ''
      if mm: mmstring = 'MULTIPLYMAPPED'
      e = self.entries[self.segments[i]['psl_entry_index']].entry
      print e['tName']+"\t"+str(e['tStart'])+"\t"+str(e['tEnd'])+"\t"+\
            e['strand']+"\t"+str(self.segments[i]['query_bed'])+"\t"+str(get_psl_quality(e))+"\t"+str(eindex)+"\t"+overstring+"\t"+mmstring
  # For the collection of alignments go through
  # all possible pairs and report any that overlap with eachother
  # in the target sequence and how much they overlap with eachother
  def find_overlapping_segment_targets(self):
    self.overlapping_segment_targets = {}
    overlapping = []
    for segindex1 in range(0,len(self.segments)):
      for segindex2 in range(segindex1+1,len(self.segments)):
        over = self.get_target_overlaps(segindex1,segindex2)
        if not over: continue
        if over[2] < self.minimum_overlap: continue #Too small to call overlapped
        overlapping.append([segindex1, segindex2, over[2]])
    for over in overlapping:
      self.overlapping_segment_targets[over[0]] = {}
      self.overlapping_segment_targets[over[1]] = {}
    for over in overlapping:
      self.overlapping_segment_targets[over[0]][over[1]] = over[2]
      self.overlapping_segment_targets[over[1]][over[0]] = over[2]
    return

  def get_target_overlaps(self,segindex1,segindex2):
    over = []
    i = self.segments[segindex1]['psl_entry_index']
    ibed = self.segments[segindex1]['query_bed']
    j = self.segments[segindex2]['psl_entry_index']
    jbed = self.segments[segindex2]['query_bed']
    ei = self.entries[i].entry
    ej = self.entries[j].entry
    iobs = set()
    for iexon in range(0,len(ei['blockSizes'])):
      for ibase in range(0,ei['blockSizes'][iexon]):
        qactual = ei['qStarts_actual'][iexon]+ibase
        t = ei['tStarts'][iexon]+ibase
        if qactual >= ibed[0] and qactual < ibed[1]:
          iobs.add(ei['tName']+':'+str(t))
    jobs = set()
    for jexon in range(0,len(ej['blockSizes'])):
      for jbase in range(0,ej['blockSizes'][jexon]):
        qactual = ej['qStarts_actual'][jexon]+jbase
        t = ej['tStarts'][jexon]+jbase
        if qactual >= jbed[0] and qactual < jbed[1]:
          jobs.add(ej['tName']+':'+str(t))
    overset = set()
    for jcoord in jobs:
      if jcoord in iobs:
        overset.add(jcoord)
    if len(overset) > 0:
      return [len(iobs),len(jobs),len(overset)]
    return False

class GenericOrderedMultipleAlignmentPSLReader():
  def __init__(self,fh=None):
    self.fh = fh
    self.previous = None
  def set_handle(self,input_fh):
    self.fh = input_fh
  def open_file(self,filename):
    self.fh = open(filename)
  def close(self):
    self.fh.close()
  def read_next(self):
    mpa = MultiplePSLAlignments()
    mcnt = 0
    current_name = None
    if self.previous:      #We have one waiting to go into an alignment
      l1 = self.previous
      p1 = PSL(l1.rstrip())
      current_name = p1.value('qName')
      mpa.add_entry(p1)
      mcnt +=  1
    else: # It must be our first entry, so prime our buffer
      l1 = None
      while True:
        l1 = self.fh.readline()
        if not l1:
          return None
        if not is_valid(l1.rstrip()): continue # go till we get a PSL
        break
      p1 = PSL(l1.rstrip())
      current_name = p1.value('qName')
      mpa.add_entry(p1)
      mcnt += 1
    while True:
      l2 = self.fh.readline()
      if not l2: 
        self.previous = None
        if mcnt > 0:
          return mpa
        return None
      if not is_valid(l2): 
        sys.stderr.write("Warning line is not a valid psl line\n"+l2.rstrip()+"\n")
        continue # just skip strange bad lines like we never saw them
      p2 = PSL(l2.rstrip())
      if p2.value('qName') == current_name: # We are working on this set of entries
        mpa.add_entry(p2)
        mcnt += 1
      else: # We have a new set so buffer it and output what we have so far
        self.previous = l2 # buffer the line
        if mcnt > 0:
          return mpa
        sys.stderr.write("ERROR: How are we here?\n")
        sys.exit()
      
def is_num(val):
  if re.match('^\d+$',str(val)): return True
  return False

# Pre: an array of PSL entries ordered by the actual query
#      So a positive strand is ready to go
#      but a negative strand set needs to be traversed backwards
#      All entries must be on the same strand and must be on the same chromosome
#         This will throw an error if not satisfied.
#      Multiple query names won't throw an error, but only the first will be used
def stitch_query_trimmed_psl_entries(entries):
  if len(entries) == 0:
    sys.stderr.write("WARNING tried stitch together zero sequences")
    return None
  strand = entries[0].value('strand')
  chrom = entries[0].value('tName')
  for e in entries:
    if e.value('strand') != strand:
      sys.stderr.write("ERROR: stitching requires same strand for all PSL")
      sys.exit()
    if e.value('tName') != chrom:
      sys.stderr.write("ERROR: stitching requires same ref sequence for all PSL")
      sys.exit()
  eordered = entries[:]
  if strand == '-':
    eordered = entries[::-1]
  prevend = 0
  outpsl = eordered[0].copy()
  tstarts = []
  qstarts = []
  bsizes = []
  for i in range(0,len(eordered)):
      #left trim by the right most value of the previous
      if eordered[i].value('tEnd') < prevend:
        sys.stderr.write("WARNING: block skipped because of order\n")
        continue
      te = eordered[i].left_t_trim(prevend+1)
      if len(tstarts) == 0:
        for j in range(0,te.value('blockCount')):
          tstarts.append(te.value('tStarts')[j])
          qstarts.append(te.value('qStarts')[j])
          bsizes.append(te.value('blockSizes')[j])
      elif tstarts[-1]+bsizes[-1]+1==te.value('tStarts')[0] and \
           qstarts[-1]+bsizes[-1]+1==te.value('qStarts')[0]:
        #Handle the special case where the next block is exactly after the previous... the are combined
        sys.stderr.write("Warning: APPEND CASE.. not a bad thing... just not common\n")
        bsizes[-1]+=te.value('blockSizes')[0]
        # The rest can be done normally
        if te.value('blockCount') > 1:
          for j in range(1,te.value('blockCount')):
            tstarts.append(te.value('tStarts')[j])
            qstarts.append(te.value('qStarts')[j])
            bsizes.append(te.value('blockSizes')[j])
      else:
        # Most normally we would just add the blocks
        for j in range(0,te.value('blockCount')):
          tstarts.append(te.value('tStarts')[j])
          qstarts.append(te.value('qStarts')[j])
          bsizes.append(te.value('blockSizes')[j]) 
      prevend = te.value('tEnd')
  outpsl.update_alignment_details(bsizes,qstarts,tstarts)
  #print len(qstarts)
  #print len(tstarts)
  #print len(bsizes)
  #print outpsl.value('blockCount')
  #print "positive strand"
  #print outpsl.get_line()
  return outpsl

