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
from SequenceBasics import rc

# A general class for a single PSL alignment
# This class will be able to perform coordinate conversions
# Pre: optionally can create with a psl line
#      alternatively you could set the entry directly
class PSL:
  def __init__(self,in_line=None):
    self.entry = None #The PSL entry in hash format when set
    self.min_intron_length = 68
    self.query_seq = None # This can optionally be set later
    self.reference_hash = None # This can optionally be set later
    if in_line:
      self.entry = line_to_entry(in_line.rstrip())
  def value(self,field_name):
    if field_name not in self.entry:
      sys.stderr.write("ERROR: "+field_name+" is not a valid field name\n")
      sys.exit()
    return self.entry[field_name]
  def set_entry(self,in_entry):
    self.entry = in_entry
  def get_entry(self):
    return self.entry
  def get_line(self):
    return entry_to_line(self.entry)
  def get_coverage(self):
    return sum(self.entry['blockSizes'])

  def pretty_print(self,window=50):
    if not self.query_seq or not self.reference_hash:
      sys.stderr.write("ERROR: Cannot pretty print unless query sequence and reference_hash have been set\n")
    exons = self.get_alignment_strings()
    z = 0
    for aligns in exons:
      z+=1
      [qexon,rexon] = aligns
      print "exon "+str(z)+'/'+str(len(exons))
      for j in range(0,len(qexon),window):
        print qexon[j:j+window]
        print rexon[j:j+window]
        print ''
  # Pre: Must already have set query_seq and reference_hash
  # Post: Returns an array [query_alignments,target_alignments]
  #       These alignments are broken apart by exon according the min_intron_length
  def get_alignment_strings(self):
    if not self.query_seq or not self.reference_hash:
      sys.stderr.write("ERROR: Cannot pretty print unless query sequence and reference_hash have been set\n")
    query = self.query_seq
    if self.value('strand') == '-': 
      query = rc(self.query_seq)
    qseq = ''
    rseq = ''
    qcur = self.value('qStarts')[0]
    rcur = self.value('tStarts')[0]
    qprev = 0
    rprev = 0
    for i in range(0,self.value('blockCount')):
      blen = self.value('blockSizes')[i]
      qS = self.value('qStarts')[i]
      qE = qS + blen
      tS = self.value('tStarts')[i]
      tE = tS + blen
      if qprev > 0 and rprev > 0: #not our first time
        qgap = qS-qprev
        tgap = tS-rprev
        #if qgap > 0 and tgap > 0:
        #  sys.stderr.write("Warning portion of psl with no alignment\n")
        if tgap-qgap >= self.min_intron_length:
          qseq += '.'
          rseq += '.'
        else:
          if qgap > 0:
            qseq += query[qE:qE+qgap].upper()
          if tgap > 0:
            rseq += self.reference_hash[self.value('tName')][tE:tE+tgap].upper()
          if qgap < tgap:
            qseq += '-'*(tgap-qgap)
          if tgap < qgap:
            rseq += '-'*(qgap-tgap)
      qseq += query[qS:qE].upper()
      rseq += self.reference_hash[self.value('tName')][tS:tE].upper()
      qprev = qE
      rprev = tE
    qexons = qseq.split('.')
    rexons = rseq.split('.')
    exons = []
    for i in range(0,len(qexons)):
      exons.append([qexons[i],rexons[i]])
    return exons

  # Pre: in_query_seq is the query sequence
  def set_query(self,in_query_seq):
    self.query_seq = in_query_seq

  # Pre: in_ref_hash is a hash for reference sequences keyed
  #      by sequence name
  def set_reference_dictionary(self,in_ref_hash):
    self.reference_hash = in_ref_hash

  # Calculate quality based on the number of mismatched bases plus the number of insertions plus the number of deletions divided by the number of aligned bases
  def get_quality(self):
    return 1-float(int(self.entry['misMatches'])+int(self.entry['qNumInsert'])+int(self.entry['tNumInsert']))/float(self.get_coverage())

  def validate(self):
    return is_valid(self.get_line())
  #Pre: 1-index coordinate
  #Post: 1-index coordinate
  def convert_coordinate_actual_query_to_target(self,aq_coord):
    qcoord = aq_coord
    if self.entry['strand'] == '-':
      qcoord = self.entry['qSize']-aq_coord+1
    for i in range(0,len(self.entry['blockSizes'])):
      for j in range(0,self.entry['blockSizes'][i]):
        if qcoord == j+self.entry['qStarts'][i]+1:
          return self.entry['tStarts'][i]+j+1
    return None
  #Pre: 1-index coordinate
  #Post: 1-index coordinate
  def convert_coordinate_target_to_actual_query(self,t_coord):
    for i in range(0,len(self.entry['blockSizes'])):
      for j in range(0,self.entry['blockSizes'][i]):
        if t_coord == j+self.entry['tStarts'][i]+1:
          if self.entry['strand'] == '-':
            return self.entry['qSize']-(self.entry['qStarts'][i]+j+1)+1
          return self.entry['qStarts'][i]+j+1
    return None
  #Pre: 1-index coordinate
  #Post: 1-index coordinate
  def convert_coordinate_query_to_actual_query(self,q_coord):
    if self.value('strand') == '+':
      return q_coord
    # must be negative strand
    return self.value('qSize')-q_coord+1

  # Pre: Have already set a psl
  # Post: Get a copy PSL object
  def copy(self):
    n = PSL()
    n.entry = self.entry.copy()
    n.min_intron_length = self.min_intron_length
    return n

  ## Remove any alignment less than the 1-index query coordiante
  #  if the strand is positive we will be cutting away the 
  #  left side because its qactual, but if its negative we will right trim
  def left_qactual_trim(self,coord):
    if self.value('strand') == '+':
      return self.left_q_trim(coord)
    return self.right_q_trim(self.value('qSize')-coord+1)
  def right_qactual_trim(self,coord):
    if self.value('strand') == '+':
      return self.right_q_trim(coord)
    return self.left_q_trim(self.value('qSize')-coord+1)
  # for these actual trimming functions trim regardess of strand and
  # return a new PSL entry.  
  # For now we lose some information on matches and mismatches
  def left_q_trim(self,incoord):
    outp = self.copy()
    qstarts = []
    bsizes = []
    tstarts = []
    for i in range(0,self.value('blockCount')):
      for j in range(0,self.value('blockSizes')[i]):
        coord0 = self.value('qStarts')[i]+j
        targ0 = self.value('tStarts')[i]+j
        if coord0+1 >= incoord: # Keep everything from this point forward
          qstarts.append(coord0) 
          tstarts.append(targ0)         
          bsizes.append(self.value('blockSizes')[i]-(coord0-self.value('qStarts')[i]))
          break # We have the start and size from this block
        #Else we are disregarding the block
    outp.update_alignment_details(bsizes,qstarts,tstarts)
    return outp
  def right_q_trim(self,incoord):
    outp = self.copy()
    qstarts = []
    bsizes = []
    tstarts = []
    for i in range(0,self.value('blockCount')):
      for j in range(0,self.value('blockSizes')[i]):
        coord0 = self.value('qStarts')[i]+j
        targ0 = self.value('tStarts')[i]+j
        if coord0+1 > incoord:
          # We have gone too far, we need to output everything we had
          # before we reached this point
          if j==0:
            outp.update_alignment_details(bsizes,qstarts,tstarts)
            return outp
          # We do have something in this loop, but recall its the previous
          bsizes.append(j)
          qstarts.append(self.value('qStarts')[i])
          tstarts.append(self.value('tStarts')[i])
          outp.update_alignment_details(bsizes,qstarts,tstarts)
          return outp
      # If we have not returned, then we can add these things to outp
      bsizes.append(self.value('blockSizes')[i])
      qstarts.append(self.value('qStarts')[i])
      tstarts.append(self.value('tStarts')[i])
    outp.update_alignment_details(bsizes,qstarts,tstarts)
    return outp

  # for these target sequence trimming functions trim based on a 1-indexed coordinate
  # for left trim delete any sequences less than this 1-indexed coordiante
  # return a new PSL entry.  
  # For now we lose some information on matches and mismatches
  def left_t_trim(self,incoord):
    outp = self.copy()
    qstarts = []
    bsizes = []
    tstarts = []
    for i in range(0,self.value('blockCount')):
      for j in range(0,self.value('blockSizes')[i]):
        coord0 = self.value('qStarts')[i]+j
        targ0 = self.value('tStarts')[i]+j
        if targ0+1 >= incoord: # Keep everything from this point forward
          qstarts.append(coord0) 
          tstarts.append(targ0)         
          bsizes.append(self.value('blockSizes')[i]-(coord0-self.value('qStarts')[i]))
          break # We have the start and size from this block
        #Else we are disregarding the block
    outp.update_alignment_details(bsizes,qstarts,tstarts)
    return outp
  def right_t_trim(self,incoord):
    outp = self.copy()
    qstarts = []
    bsizes = []
    tstarts = []
    for i in range(0,self.value('blockCount')):
      for j in range(0,self.value('blockSizes')[i]):
        coord0 = self.value('qStarts')[i]+j
        targ0 = self.value('tStarts')[i]+j
        if targ0+1 > incoord:
          # We have gone too far, we need to output everything we had
          # before we reached this point
          if j==0:
            outp.update_alignment_details(bsizes,qstarts,tstarts)
            return outp
          # We do have something in this loop, but recall its the previous
          bsizes.append(j)
          qstarts.append(self.value('qStarts')[i])
          tstarts.append(self.value('tStarts')[i])
          outp.update_alignment_details(bsizes,qstarts,tstarts)
          return outp
      # If we have not returned, then we can add these things to outp
      bsizes.append(self.value('blockSizes')[i])
      qstarts.append(self.value('qStarts')[i])
      tstarts.append(self.value('tStarts')[i])
    outp.update_alignment_details(bsizes,qstarts,tstarts)
    return outp

  def update_alignment_details(self,blocksizes,qstarts,tstarts):
    #Take in new alignment details and update the psl
    #Recalculating the stats will wipe some mismatch details
    self.entry['qStarts'] = qstarts
    self.entry['tStarts'] = tstarts
    self.entry['blockSizes'] = blocksizes
    self.entry['qStart'] = qstarts[0]
    self.entry['qEnd'] = qstarts[-1]+blocksizes[-1]
    self.entry['tStart'] = tstarts[0]
    self.entry['tEnd'] = tstarts[-1]+blocksizes[-1]
    self.entry['blockCount'] = len(blocksizes)
    self.entry['qStarts_actual'] = calculate_qstarts_actual(self.value('qSize'),qstarts,blocksizes,self.value('strand'))
    self.recalculate_stats()

  #Pre: A psl entry
  #Post: Forget any infromation about the actual sequence
  #      and recalculate stats
  def recalculate_stats(self):
    self.entry['matches'] = 0
    self.entry['misMatches'] = 0
    self.entry['repMatches'] = 0
    self.entry['nCount'] = 0
    self.entry['qNumInsert'] = 0
    self.entry['qBaseInsert'] = 0
    self.entry['tNumInsert'] = 0
    self.entry['tBaseInsert'] = 0
    for i in range(0,len(self.entry['blockSizes'])):
      for j in range(0,self.entry['blockSizes'][i]):
        self.entry['matches'] += 1
    if len(self.entry['blockSizes']) < 2: return
    for i in range(1,len(self.entry['blockSizes'])):
      qcur = self.entry['qStarts'][i]
      qpast = self.entry['qStarts'][i-1]+self.entry['blockSizes'][i-1]
      qgap = qcur-qpast
      if qgap > 0:
        self.entry['qNumInsert']+=1
        self.entry['qBaseInsert']+=qgap
      tcur = self.entry['tStarts'][i]
      tpast = self.entry['tStarts'][i-1]+self.entry['blockSizes'][i-1]
      tgap = tcur-tpast
      if tgap > 0 and tgap < self.min_intron_length:
        self.entry['tNumInsert']+=1
        self.entry['tBaseInsert']+=tgap
    return

  #Pre: A psl entry
  #     self.query_seq and self.reference_hash must be set
  #Post: Use any infromation about the actual sequence
  #      and recalculate stats
  def correct_stats(self):
    if not self.query_seq or not self.reference_hash:
      sys.stderr.write("ERROR: Cannot recalculate stats on the psl without both the query and reference sequences\n")
      sys.exit()
    g = self.reference_hash
    query = self.query_seq
    if self.value('strand') == '-':
      query = rc(self.query_seq)
    nCount = 0
    matches = 0
    misMatches = 0
    prev_qE = 0
    prev_tE = 0
    qNumInsert = 0
    qBaseInsert = 0
    tNumInsert = 0
    tBaseInsert = 0
    for i in range(self.value('blockCount')):
      blen = self.value('blockSizes')[i]
      qS = self.value('qStarts')[i] #query start
      qE = qS + blen             #query end
      tS = self.value('tStarts')[i] #target start
      tE = tS + blen             #target end
      #Work on gaps
      if prev_qE > 0 or prev_tE > 0: #if its not our first time through
        tgap = tS-prev_tE
        if tgap < self.min_intron_length and tgap > 0:
          tNumInsert += 1
          tBaseInsert += tgap
        qgap = qS-prev_qE
        if qgap > 0:
          qNumInsert += 1
          qBaseInsert += qgap
      qseq = query[qS:qE].upper()
      rseq = g[self.value('tName')][tS:tE].upper()
      #print qseq+"\n"+rseq+"\n"
      for j in range(0,blen):
        if qseq[j] == 'N':
          nCount += 1
        elif qseq[j] == rseq[j]:
          matches += 1
        else:
          misMatches += 1
      prev_qE = qE
      prev_tE = tE
    self.entry['matches'] = matches
    self.entry['misMatches'] = misMatches
    self.entry['nCount'] = nCount
    self.entry['qNumInsert'] = qNumInsert
    self.entry['qBaseInsert'] = qBaseInsert
    self.entry['tNumInsert'] = tNumInsert
    self.entry['tBaseInsert'] = tBaseInsert
    self.entry['qSize'] = len(query)
    self.entry['tSize'] = len(g[self.value('tName')])

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
  #return 1-float(int(e['misMatches'])+int(e['qNumInsert'])+int(e['tNumInsert']))/float(get_coverage(e))
  return 1-float(int(e['misMatches'])+int(e['qBaseInsert'])+int(e['tBaseInsert']))/float(get_coverage(e))


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
  v['matches'] = int(f[0])
  v['misMatches'] = int(f[1])
  v['repMatches'] = int(f[2])
  v['nCount'] = int(f[3])
  v['qNumInsert'] = int(f[4])
  v['qBaseInsert'] = int(f[5])
  v['tNumInsert'] = int(f[6])
  v['tBaseInsert'] = int(f[7])
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
  v['qStarts_actual'] = calculate_qstarts_actual(v['qSize'],v['qStarts'],v['blockSizes'],v['strand'])
  return v

def calculate_qstarts_actual(qSize,qStarts,blockSizes,strand):
  if strand == '+':
    return qStarts # making life easier
  qStarts_actual = []
  for i in range(0,len(blockSizes)):
    qStarts_actual.append(qSize-(qStarts[i]+blockSizes[i]))    
  return qStarts_actual

class MultiplePSLAlignments:
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
  def __init__(self):
    self.fh = None
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
      
def is_valid(line):
  # Test if a PSL line is valid.
  f = line.rstrip().split("\t")
  if len(f) > 21:
    sys.stderr.write("Error: Line is longer than 21 fields\n")
    return False
  if len(f) < 21:
    sys.stderr.write("Error: Line is less than 21 fields\n")
    return False
  if not is_num(f[0]):
    sys.stderr.write("Error: matches is not numeric\n")
    return False
  if not is_num(f[1]):
    sys.stderr.write("Error: misMatches is not numeric\n")
    return False
  if not is_num(f[2]):
    sys.stderr.write("Error: repMatches is not numeric\n")
    return False
  if not is_num(f[3]):
    sys.stderr.write("Error: nCount is not numeric\n")
    return False
  if not is_num(f[4]):
    sys.stderr.write("Error: qNumInsert is not numeric\n")
    return False
  if not is_num(f[5]):
    sys.stderr.write("Error: qBaseInsert is not numeric\n")
    return False
  if not is_num(f[6]):
    sys.stderr.write("Error: tNumInsert is not numeric\n")
    return False
  if not is_num(f[7]):
    sys.stderr.write("Error: tBaseInsert is not numeric\n")
    return False
  if not re.match('^[+-]$',f[8]):
    sys.stderr.write("Error: strand not + or -\n")
    return False
  if not is_num(f[10]):
    sys.stderr.write("Error: qSize is not numeric\n")
    return False
  e = line_to_entry(line)
  if len(e['blockSizes']) != len(e['qStarts']) or len(e['blockSizes']) != len(e['tStarts']):
    sys.stderr.write("Error: Block sizes is not the same as query or target size\n")
    return False
  if len(e['qStarts']) > 1:
    for i in range(1,len(e['qStarts'])):
      if e['qStarts'][i] <= e['qStarts'][i-1]:
        sys.stderr.write("Error: Query starts are not in ascending order\n")  
        return False
  if len(e['tStarts']) > 1:
    for i in range(1,len(e['tStarts'])):
      if e['tStarts'][i] <= e['tStarts'][i-1]:
        sys.stderr.write("Error: Target starts are not in ascending order\n")  
        return False
  return True

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
