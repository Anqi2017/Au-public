# psl_basics.py
#
# A module for holding functions for handling psl files
# 
# These include:
# line_to_entry() - read a line from a psl file into a dictionary
#           also includes a new array of adjusted coordinates 
#           in the same coordiante system as a positive strand
#           for easier comparison of query alignments
#
# query_base_overlap_size() - for determining how many bases of two
#           psl entries map to the same query location.  This is helpful
#           in determining if two psl entries from the same query
#           are describing different parts of the same query. For a faster
#           calculation see query_region_overlap_size().
# query_coordinates_base_overlap_size() - like above but does not assume
#           the queries should have the same query name.  For special cases
#           where the query name has been modified.
#
# query_region_overlap_size() - for determining how many bases overlap within
#           the overall region spanned accross the query between two psl
#           entries.  This is faster than the base-wise calculation.
# query_coordinates_base_overlap_size() - like above but does not assume 
#           the queries should have the same query name.  For special cawses
#           where the query name has been modified.
#
# query_gap_size() - for determining the size of the gap between two psl
#           entries on the same query.  
# query_coordinates_gap_size() - like above but does not assume
#           the queries should have the same query name.  For special cawses
#           where the query name has been modified.
#
# convert_entry_to_genepred_line() - change an entry (in dictionary format)
#           to a genepred format line
import re, sys
from RangeBasics import GenomicRange

# pre: one psl entry, one coordinate (1-indexed)
# pos: one coordinate (1-indexed)
def convert_target_to_query_coord(psl,coord):
  for i in range(0,len(psl['blockSizes'])):
    z = 0
    for j in range(psl['tStarts'][i],psl['tStarts'][i]+psl['blockSizes'][i]):
      if j == coord-1:
        return psl['qStarts_actual'][i]+z+1
      z += 1
  return

def convert_entry_to_target_bed(psl,color):
  ostring = psl['tName'] + "\t" #chrom
  ostring += str(psl['tStart']) + "\t" #chromStart
  ostring += str(psl['tEnd']) + "\t" #chromEnd
  ostring += psl['qName'] + "\t" #name
  ostring += "1" + "\t" #score
  ostring += psl['strand'] + "\t" #strand
  ostring += str(psl['tStart']) + "\t" #thickStart
  ostring += str(psl['tEnd']) + "\t" #thickEnd
  ostring += color + "\t" #itemRgb
  ostring += str(psl['blockCount']) + "\t" #blockCount
  ostring += ",".join([str(x) for x in psl['blockSizes']])+"," + "\t"
  ostring += ",".join([str(x-psl['tStart']) for x in psl['tStarts']])+","
  return ostring

# pre: one psl entry
# post: one genepred line of the query's location on the target
#       we will use the genepred format with an aditional gene name field
def convert_entry_to_genepred_line(psl):
  vals = [psl['qName'], psl['qName'], psl['tName'], psl['strand'], str(psl['tStart']), str(psl['tEnd']), str(psl['tStart']), str(psl['tEnd']), str(psl['blockCount'])]
  starts = []
  ends = []
  if len(psl['blockSizes']) != len(psl['tStarts']):
    sys.stderr.write("ERROR different sizes and starts for \n"+str(psl)+"\n")
    sys.stderr.write(str(len(psl['blockSizes']))+"\t"+str(len(psl['tStarts']))+"\n")
    return False
  for i in range(0,len(psl['blockSizes'])):
    starts.append(str(psl['tStarts'][i]))
    ends.append(str(psl['tStarts'][i]+psl['blockSizes'][i]))
  vals.append(','.join(starts)+',')
  vals.append(','.join(ends)+',')
  return "\t".join(vals)

# pre: two psl entries 
# post: return the size of the gap (number of bases in the gap), 
#       zero if there is no gap
#         or if they occur on different query names
#         of if they overlap its still zero
def query_gap_size(e1,e2):
  if e1['qName'] != e2['qName']:
    return 0
  return query_coordinates_gap_size(e1,2)

# pre: two psl entries 
# post: return the size of the gap, 
#       zero if there is no gap, or there is overlap
#       checks query coordinates regardless of the name
# modifies: none
def query_coordinates_gap_size(e1,e2):
  if e2['qStart'] > e1['qEnd']-1:
    # they start is zero indexed and end is one indexed
    dist = e2['qStart']-(e1['qEnd']-1)
    # gap size is 1 less than the distance
    return dist - 1
  if e1['qStart'] > e2['qEnd']-1:
    dist = e1['qStart']-(e2['qEnd']-1)
    return dist - 1
  return 0

# pre: two psl entries 
# post: see if the queries overlap and return how many bases
#       zero if there is no overlap or they are different names
# modifies: none
def query_region_overlap_size(e1,e2):
  if e1['qName'] != e2['qName']:
    return 0
  return query_coordinates_region_overlap_size(e1,2)

# pre: two psl entries 
# post: see if the queries overlap and return how many bases
#       zero if there is no overlap or they are different names
# modifies: none
def query_base_overlap_size(e1,e2):
  if e1['qName'] != e2['qName']:
    return 0
  return query_coordinates_base_overlap_size(e1,2)

# pre: two psl entries
# post: see if the query coordinates overlap for the two entries and count how many bases overlap
#       we don't check the names of queries here for special cases where you want to compare them regardless of what has happened to the query name
#       zero if there is no overlap, otherwise number of bases that overlap
# modifies: none
def query_coordinates_region_overlap_size(e1,e2):
  e1start = e1['qStart']
  e1end = e1['qEnd']-1
  e2start = e2['qStart']
  e2end = e2['qEnd']-1
  # not matching
  if e1end < e2start:
    return 0
  if e2end < e1start:
    return 0
  #encompassing
  if e1start <= e2start and e1end >= e2end:
    return e2end - e2start + 1
  if e2start <= e1start and e2end >= e1end:
    return e1end - e1start + 1
  #overlapping left
  if e1start < e2start and e1end >= e2start:
    return e1end - e2start + 1
  if e2start < e1start and e2end >= e1start:
    return e2end - e1start + 1
  #overlapping right
  if e1start <= e2end and e1end > e2end:
    return e2end - e1start + 1
  if e2start <= e1end and e2end > e1end:
    return e1end - e2start + 1
  print 'ERROR psl_basics'
  sys.exit()

# pre: two psl entries
# post: see if the query coordinates overlap for the two entries and count how many bases overlap
#       we don't check the names of queries here for special cases where you want to compare them regardless of what has happened to the query name
#       zero if there is no overlap, otherwise number of bases that overlap
# modifies: none
def query_coordinates_base_overlap_size(e1,e2):
  # to speed things up if there is no region overlap, then don't look for a base overlap
  if query_coordinates_region_overlap_size(e1,e2) == 0: return 0
  e1bases = set()
  for i in range(0,len(e1['blockSizes'])):
    for j in range(e1['qStarts_actual'][i],e1['qStarts_actual'][i]+e1['blockSizes'][i]):
      e1bases.add(j)
  baseoverlap = 0
  for i in range(0,len(e2['blockSizes'])):
    for j in range(e2['qStarts_actual'][i],e2['qStarts_actual'][i]+e2['blockSizes'][i]):
      if j in e1bases:
        baseoverlap += 1
  return baseoverlap

# Put our entry stored PSL back into a PSL line
def entry_to_line(e):
  line = str(e['matches']) + "\t" + str(e['misMatches']) + "\t" \
       + str(e['repMatches']) + "\t" + str(e['nCount']) + "\t" \
       + str(e['qNumInsert']) + "\t" + str(e['qBaseInsert']) + "\t" \
       + str(e['tNumInsert']) + "\t" + str(e['tBaseInsert']) + "\t" \
       + e['strand'] + "\t" + e['qName'] + "\t" + str(e['qSize']) + "\t" \
       + str(e['qStart']) + "\t" + str(e['qEnd']) + "\t" + e['tName'] + "\t" \
       + str(e['tSize']) + "\t" + str(e['tStart']) + "\t" \
       + str(e['tEnd']) + "\t" + str(e['blockCount']) + "\t" \
       + ','.join([str(x) for x in e['blockSizes']])+ ',' + "\t" \
       + ','.join([str(x) for x in e['qStarts']])+',' + "\t" \
       + ','.join([str(x) for x in e['tStarts']])+','
  return line

def get_coverage(e):
  return sum(e['blockSizes'])

# Calculate quality based on the number of mismatched bases plus the number of insertions plus the number of deletions divided by the number of aligned bases
def get_quality(e):
  return 1-float(int(e['misMatches'])+int(e['qNumInsert'])+int(e['tNumInsert']))/float(get_coverage(e))


# pre: a line from a psl file
# post: a dictionary with terms assigned for each variable in the psl line
#       includes an additional 'qStarts_actual' keyed array that contains the
#         start locations for the query in the same coordinate system as the
#         it would have on the positive strand to simplify comparisions of 
#         coordinate overlaps between multiple hits of the same query regardless
#         of direction.
# modifies: none
def line_to_entry(line):
  line = line.strip()
  f = line.split("\t")
  v = {}
  v['matches'] = f[0]
  v['misMatches'] = f[1]
  v['repMatches'] = f[2]
  v['nCount'] = f[3]
  v['qNumInsert'] = f[4]
  v['qBaseInsert'] = f[5]
  v['tNumInsert'] = f[6]
  v['tBaseInsert'] = f[7]
  v['strand'] = f[8]
  if not re.match('^[+-]$',v['strand']):
    print "unsupported strand type "+v['strand']
    sys.exit()
  v['qName'] = f[9]
  v['qSize'] = int(f[10])
  v['qStart'] = int(f[11])
  v['qEnd'] = int(f[12])
  v['tName'] = f[13]
  v['tSize'] = int(f[14])
  v['tStart'] = int(f[15])
  v['tEnd'] = int(f[16])
  v['blockCount'] = int(f[17])
  v['blockSizes'] = map(int,f[18].strip(',').split(','))
  v['qStarts'] = map(int,f[19].strip(',').split(','))
  v['tStarts'] = map(int,f[20].strip(',').split(','))
  v['qStarts_actual'] = v['qStarts'] # making life easier
  if v['strand'] == '-':
    v['qStarts_actual'] = []
    for i in range(0,len(v['blockSizes'])):
      v['qStarts_actual'].append(v['qSize']-(v['qStarts'][i]+v['blockSizes'][i]))    
  return v

class MultiplePSLAlignments:
  def __init__(self):
    self.entries = []
    self.minimum_coverage = 1 #how many base pairs an alignment must cover to be part of a multiple alignment
    self.qName = None
    self.best_coverage_fraction = 0.9 #how much of an alignment be the best alignment
                                       #where it aligns to be kept later on
    self.multiply_mapped_minimum_quality_difference = 0.05 # if the qualities are similar we can
                                                           # treat it as an ambiguous alignment
    self.multiply_mapped_minimum_overlap = 50 # The minimum number of bases that must be shared between two alignments to consider them multiply mapped
    self.best_alignment = None # This will hold a BestAlignmentCollection class and is set by 'best_query'
    self.verbose = False
  def set_minimum_coverage(self,mincov):
    self.minimum_coverage = mincov
  def add_line(self,line):
    self.add_entry(line_to_entry(line))
  def add_entry(self,entry):
    if not self.qName:
      self.qName = entry['qName']
    else:
      if entry['qName'] != self.qName:
        sys.stderr.write("WARNING multiple alignments must have the same query name.  This entry will not be added\n")
        return False
    if self.minimum_coverage > 1:
      cov = get_coverage(entry)
      if cov < self.minimum_coverage:
        if self.verbose: sys.stderr.write("WARNING alignment less than minimum coverage.\n")
        return False
    self.entries.append(entry)
    return True
  def get_alignment_count(self):
    return len(self.entries)
  def get_tNames(self):
    names = set()
    for name in [x['tName'] for x in self.entries]:
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
      e = self.entries[eindex]
      # Go through each block of the alignment
      for i in range(0,e['blockCount']):
        # Set relevant mapped alignments for each base of the query
        for z in range(e['qStarts_actual'][i],e['qStarts_actual'][i]+e['blockSizes'][i]):
          if z not in query: query[z] = {}
          query[z][eindex] = {}
          query[z][eindex]['tName'] = e['tName']
          query[z][eindex]['coverage'] = get_coverage(e)
          query[z][eindex]['quality'] = get_quality(e)
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
      es[i]['coverage'] = get_coverage(self.entries[i])      
      es[i]['quality'] = get_quality(self.entries[i])      
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

# Store the result of a 'best_query' this
class BestAlignmentCollection:
  def __init__(self):
    self.entries = {}  # psl entries stored by an integer key
    self.segments = [] # contains a query_bed and a psl_entry_index
    self.qName = None
    self.minimum_overlap = 1 # by default consider any overlap as reportable overlap
    self.overlapping_segment_targets = None # set by find_overlapping_segment_targets
    self.minimum_locus_distance = 400000 # minimum number of bases to consider something a different locus 
    return

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

  def locus_count(self):
    loci = []
    for i in range(0,len(self.segments)):
      loci.append(set([i]))
    prev_count = -1
    while len(loci) != prev_count:
      #Try to combine down loci
      prev_count = len(loci)
      loci = self.combine_down(loci)
    return len(loci)

  def combine_down(self,loci):
    new_loci = []
    for i in range(0,len(loci)):
      nset = set()
      chr1 = ''
      s1 = 100000000000
      f1 = 0
      for l in loci[i]:  
        el = self.entries[self.segments[l]['psl_entry_index']]
        nset.add(l)
        chr1 = el['tName']
        if el['tStart'] < s1: s1 = el['tStart']
        if el['tEnd'] > f1: f1 = el['tEnd']
      new_loci.append(nset)
      for j in range(i+1,len(loci)):
        chr2 = ''
        s2 = 100000000000
        f2 = 0
        mset = set()
        for m in loci[j]:
          em = self.entries[self.segments[j]['psl_entry_index']]
          mset.add(m)
          chr2 = em['tName']
          if em['tStart'] < s2: s2 = em['tStart']
          if em['tEnd'] > f2: f2 = em['tEnd'] 
        gri = GenomicRange(chr1,s1+1,f1)
        grj = GenomicRange(chr2,s2+1-self.minimum_locus_distance,f2+self.minimum_locus_distance)
        if gri.overlaps(grj):
          for z in mset:
            new_loci[-1].add(z)
          for k in range(j+1,len(loci)):
            new_loci.append(loci[k])
          return new_loci
    return new_loci

  def print_report(self):
    if not self.overlapping_segment_targets:
      self.find_overlapping_segment_targets()
    print '-----'
    print self.qName
    print str(self.locus_count())+" loci"
    for i in range(0,len(self.segments)):
      overstring = ''
      if i in self.overlapping_segment_targets: overstring = 'OVERLAPPED'
      eindex = self.segments[i]['psl_entry_index']
      mm = self.segments[i]['multiply_mapped']
      e = self.entries[self.segments[i]['psl_entry_index']]
      print e['tName']+"\t"+str(e['tStart'])+"\t"+str(e['tEnd'])+"\t"+\
            e['strand']+"\t"+str(self.segments[i]['query_bed'])+"\t"+str(eindex)+"\t"+overstring+"\t"+str(mm)
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
    ei = self.entries[i]
    ej = self.entries[j]
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

