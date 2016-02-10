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
from RangeBasics import GenomicRange, Bed
from SequenceBasics import rc as rc_seq

# A general class for a single PSL alignment
# This class will be able to perform coordinate conversions
# Pre: optionally can create with a psl line
#      alternatively you could set the entry directly
class PSL:
  def __init__(self,in_line=None):
    self.entry = None #The PSL entry in hash format when set
    self.min_intron_length = 68
    self.query_seq = None # This can optionally be set later
    self.quality_seq = None
    self.reference_hash = None # This can optionally be set later
    if in_line:
      self.entry = line_to_entry(in_line.rstrip())
  def get_genepred_line(self):
    return convert_entry_to_genepred_line(self.entry)
  def value(self,field_name):
    if field_name not in self.entry:
      sys.stderr.write("ERROR: "+field_name+" is not a valid field name\n")
      sys.exit()
    return self.entry[field_name]

  def get_target_bed(self):
    return Bed(self.value('tName'),self.value('tStart'),self.value('tEnd'),self.value('strand'))

  # Returns the actual query coverage
  def get_query_bed(self):
    s1 = self.value('qStarts_actual')[0]
    s2 = self.value('qStarts_actual')[-1]+self.value('blockSizes')[-1]
    if self.value('strand') == '-':
      s1 = self.convert_coordinate_query_to_actual_query(self.value('qStarts')[-1]+self.value('blockSizes')[-1])-1
      s2 = self.convert_coordinate_query_to_actual_query(self.value('qStarts')[0]+1)
    return Bed(self.value('qName'),s1,s2)


  # Pre: Another PSL entry, whether or not to use the direction
  # Post: Return the distance between them or 0 if overlapping, or -1 if same chromosome
  def target_distance(self,psl_entry,use_direction=False):
    if self.value('tName') != psl_entry.value('tName'):
      return -1
    if use_direction and self.value('strand') != psl_entry.value('strand'):
      return -1
    b1 = Bed(self.entry['tName'],self.entry['tStart'],self.entry['tEnd'])
    b2 = Bed(psl_entry.entry['tName'],psl_entry.entry['tStart'],psl_entry.entry['tEnd'])
    if b1.overlaps(b2):
      return 0
    if b1.end < b2.start:
      return b2.start-b1.end-1
    if b1.start > b2.end:
      return b1.start-b2.end-1
    sys.stderr.write("ERROR un accounted for state\n")
    sys.exit()

  # Pre: Another PSL entry, optionally whether or not to use the direction
  # Post: Return the distance between them or 0 if overlapping, or -1 if not same chromosome
  def query_distance(self,psl_entry):
    if self.value('qName') != psl_entry.value('qName'):
      return -1
    #if use_strand and self.value('strand') != psl_entry.value('strand'):
    #  return -1
    b1 = self.get_query_bed()
    b2 = psl_entry.get_query_bed()
    if b1.overlaps(b2):
      return 0
    if b1.end < b2.start:
      return b2.start-b1.end-1
    if b1.start > b2.end:
      return b1.start-b2.end-1
    sys.stderr.write("ERROR un accounted for state\n")
    sys.exit()

  def target_overlap_size(self,psl2,use_direction=False):
    if self.value('tName') != psl2.value('tName'):
      return 0
    if use_direction and self.value('strand') != psl2.value('strand'):
      return 0
    # on same chromosome
    output = 0
    for i in range(0,self.value('blockCount')):
      for j in range(0,psl2.value('blockCount')):
        b1 = Bed(self.value('tName'),self.value('tStarts')[i],self.value('tStarts')[i]+self.value('blockSizes')[i])
        b2 = Bed(psl2.value('tName'),psl2.value('tStarts')[j],psl2.value('tStarts')[j]+psl2.value('blockSizes')[j])
        size = b1.overlap_size(b2)
        output += size
    return output
  def query_overlap_size(self,psl2):
    if self.value('qName') != psl2.value('qName'):
      return 0
    # on same query
    output = 0
    for i in range(0,self.value('blockCount')):
      for j in range(0,psl2.value('blockCount')):
        b1 = Bed(self.value('qName'),self.value('qStarts_actual')[i],self.value('qStarts_actual')[i]+self.value('blockSizes')[i])
        b2 = Bed(psl2.value('qName'),psl2.value('qStarts_actual')[j],psl2.value('qStarts_actual')[j]+psl2.value('blockSizes')[j])
        size = b1.overlap_size(b2)
        output += size
    return output

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
      query = rc_seq(self.query_seq)
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

  def get_query(self):
    return self.query_seq
  # Pre: in_query_seq is the query sequence
  def set_quality_seq(self,in_quality_seq):
    self.quality_seq = in_quality_seq

  def get_quality_seq(self):
    return self.quality_seq

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
    n.query_seq = self.query_seq
    n.quality_seq = self.quality_seq 
    n.reference_hash = self.reference_hash
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
      query = rc_seq(self.query_seq)
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
  # Reverse complement the query and its alignment
  # Post: return a new PSL object with the reverse complemented query alignment
  def rc(self):
    p = self.copy()
    if self.value('strand')=='+':
      p.entry['strand'] = '-'
    else:
      p.entry['strand'] = '+'
    if self.query_seq:  p.query_seq = rc_seq(self.query_seq)
    if self.quality_seq: p.quality_seq = self.quality_seq[::-1]
    return p

  def concatonate(self,p2):
    #alignments are on the same chromosome and strand
    if self.value('tName') != p2.value('tName') or self.value('strand') != p2.value('strand'):
      sys.stderr.write('ERROR: cant concatonate on different chromosomes or strands')
      sys.exit()
    out = self.copy()
    if out.query_distance(p2) == 0: # overlaps queries
      out = out.right_q_trim(p2.value('qStarts')[0]+1)
      print 'overlaps query'
    if out.target_distance(p2) == 0: #overlaps targets
      out = out.right_t_trim(p2.value('tStarts')[0]+1)
      print 'overlaps target'
    l1 = out.get_line().rstrip().split("\t")
    l2 = p2.get_line().rstrip().split("\t")
    #fix blocks
    # add to the right side
    l1[18] = l1[18]+l2[18]
    l1[19] = l1[19]+l2[19]
    l1[20] = l1[20]+l2[20]
    l1[17] = str(int(l1[17])+int(l2[17]))
    l1[16] = str(int(l2[16]))
    return PSL("\t".join(l1))

  # Pre: a PSL entries on the same chromosome and strand
  #      assume the overlap of the target could be described equally well by either mate
  #      return a new psl entry where queries have been concatonated and 
  #      a new query sequence formed if query sequences are available
  #      Used for making concatonated short read alignments
  # Post Contatonates self (left) with the input argument (right)
  def concatonate_queries(self,p2):
    #p2r = p2.rc()
    if self.value('tName') != p2.value('tName') or self.value('strand') != p2.value('strand'):
      sys.stderr.write("ERROR cant concatonate if not on same chrom and strand\n")
      sys.exit()
    if self.value('tStart') == p2.value('tStart') or self.value('tEnd') == p2.value('tEnd'):
      if self.get_coverage() > self.get_coverage():
        return self.copy()
      else: return p2.copy()
    #Check overlapping cases
    if self.value('tStart') < p2.value('tStart') and self.value('tEnd') > p2.value('tEnd'):
      return self.copy()
    if p2.value('tStart') < self.value('tStart') and p2.value('tEnd') > self.value('tEnd'):
      return p2.copy()
    if self.value('strand') == '+' and self.value('tStart') > p2.value('tStart'):
      sys.stderr.write("ERROR: Unexpected order\n")
      sys.exit()
    if self.value('strand') == '-' and self.value('tStart') < p2.value('tStart'):
      sys.stderr.write("ERROR: Unexpected order\n")
      sys.exit()
    # lets put things together
    # First lets decide on a p1 and a p2
    p1 = self.copy()
    if self.value('tStart') > p2.value('tStart'):
      p1 = p2.copy()
      p2 = self.copy()
      sys.stderr.write("Warning concatonating a sequence that is to the left\n")
    # now p1 always starts first
    # see if there is a gap betwen them
    p1trim = p1.right_t_trim(p2.value('tStart'))
    if p1.query_seq:
      p1trim.query_seq = p1.query_seq[0:p1trim.value('qEnd')]
      p1trim.quality_seq = p1.quality_seq[0:p1trim.value('qEnd')]
      p1trim.entry['qSize'] = len(p1trim.query_seq) #force downsize to the mapped trimmed part
      if p1.value('strand') == '-':
        p1trim.query_seq = rc_seq(rc_seq(p1.query_seq)[0:p1trim.value('qEnd')])
        p1trim.quality_seq = ((p1.quality_seq[::-1])[0:p1trim.value('qEnd')])[::-1]
        p1trim.entry['qSize'] = len(p1trim.query_seq) #force downsize to the mapped trimmed part
    p1 = p1trim
    output = p1.copy()
    if p1.query_seq and p2.query_seq:
      new_query = p1.query_seq + p2.query_seq
      new_quality = p1.quality_seq + p2.quality_seq
      if p1.value('strand') == '-':
        new_query = rc_seq(p1.query_seq)+rc_seq(p2.query_seq)
        new_quality = p1.quality_seq[::-1]+p2.quality_seq[::-1]
      output.set_query(new_query)
      output.set_quality_seq(new_quality)
    output.entry['qSize'] = len(new_query)
    #else:
    #  output.entry['qSize'] = p1.value('qSize')+p2.value('qSize')
    new2qstarts = [x+p1.value('qSize') for x in p2.value('qStarts')]
    # First handle the case of the first p2 entry if no gap
    if new2qstarts[0] == p1.value('qEnd') and p2.value('tStart') == p1.value('tEnd'): 
      output.entry['blockSizes'][-1] += p2.value('blockSizes')[0]
    else:
      output.entry['qStarts'].append(new2qstarts[0])
      output.entry['tStarts'].append(p2.value('tStarts')[0])
      output.entry['blockSizes'].append(p2.value('blockSizes')[0])
    if len(new2qstarts) > 1:
      for i in range(1,len(new2qstarts)):
        output.entry['qStarts'].append(new2qstarts[i])
        output.entry['tStarts'].append(p2.value('tStarts')[i])
        output.entry['blockSizes'].append(p2.value('blockSizes')[i])
    output.update_alignment_details(output.value('blockSizes'),output.value('qStarts'),output.value('tStarts'))
    return output

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


def get_psl_quality(entry):
  return float(entry['matches'])/float(entry['matches']+entry['misMatches']+entry['tNumInsert']+entry['qNumInsert'])

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

#Pre: a psl entry
#Post: query_beds and target_beds, two arrays of beds
def get_beds_from_entry(entry,use_direction=False):
  query_beds = []
  target_beds = []
  print entry
  for i in range(0,entry['blockCount']):
    if use_direction:
      tb = Bed(entry['tName'],entry['tStarts'][i],entry['tStarts'][i]+entry['blockSizes'][i],entry['strand'])
      target_beds.append(tb)
    else:
      tb = Bed(entry['tName'],entry['tStarts'][i],entry['tStarts'][i]+entry['blockSizes'][i])
      target_beds.append(tb)
    qb = Bed(entry['qName'],entry['qStarts_actual'][i],entry['qStarts_actual'][i]+entry['blockSizes'][i])
    query_beds.append(tb)
  return [query_beds, target_beds]
