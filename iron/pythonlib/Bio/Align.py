class Alignment:
  def __init__(self):
    self._alignment_ranges = None
    self._query_sequence = None
    self._query_length = 0
    self._target_length = 0
    sef._set_alignment_ranges()
    return

  # These methods need to be overridden by an alignment type
  def _set_alignment_ranges(self):
    self._alignment_ranges = None
  def get_query_sequence(self):
    return self._query_sequence
  def get_target_length(self):
    return self._target_length
  def get_strand(self):
    return self._query_direction

  # These methods may be overrriden by an alignment type to just return themself
  def get_PSL(self):
    print 'hello'
    for r in self._alignment_ranges:
      print r[0].get_range_string()
      print r[1].get_range_string()
    matches = sum([x[0].length() for x in self._alignment_ranges]) # 1. Matches - Number of matching bases that aren't repeats
    misMatches = 0 # 2. Mismatches - Number of baess that don't match
    repMatches = 0 # 3. repMatches - Number of matching baess that are part of repeats
    nCount = 0 # 4. nCount - Number of 'N' bases
    qNumInsert = 0 # 5. qNumInsert - Number of inserts in query
    qBaseInsert = 0 # 6. qBaseInsert - Number of bases inserted into query
    tNumInsert = 0 # 7. Number of inserts in target
    tBaseInsert = 0 # 8. Number of bases inserted into target
    strand = self.get_strand() # 9. strand 
    qName = self._alignment_ranges[0][1].chr # 10. qName - Query sequence name
    qSize = self.get_query_length()
    qStart = self._alignment_ranges[0][1].start-1
    qEnd = self._alignment_ranges[-1][1].end
    tName = self._alignment_ranges[0][0].chr
    tSize = self.get_target_length()
    tStart = self._alignment_ranges[0][0].start-1
    tEnd = self._alignment_ranges[-1][0].end
    blockCount = len(self._alignment_ranges)
    blockSizes = ','.join([str(x[0].length()) for x in self._alignment_ranges])+','
    qStarts = ','.join([str(x[1].start-1) for x in self._alignment_ranges])+','
    tStarts = ','.join([str(x[0].start-1) for x in self._alignment_ranges])+','

    psl_string = str(matches)+"\t"+\
    str(misMatches)+"\t"+\
    str(repMatches)+"\t"+\
    str(nCount)+"\t"+\
    str(qNumInsert)+"\t"+\
    str(qBaseInsert)+"\t"+\
    str(tNumInsert)+"\t"+\
    str(tBaseInsert)+"\t"+\
    strand+"\t"+\
    qName+"\t"+\
    str(qSize)+"\t"+\
    str(qStart)+"\t"+\
    str(qEnd)+"\t"+\
    tName+"\t"+\
    str(tSize)+"\t"+\
    str(tStart)+"\t"+\
    str(tEnd)+"\t"+\
    str(blockCount)+"\t"+\
    blockSizes+"\t"+\
    qStarts+"\t"+\
    tStarts
    print psl_string
    return

  def get_SAM(self):
    return
