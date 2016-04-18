from Bio.Sequence import rc
from string import maketrans
class Alignment:
  def __init__(self):
    self._alignment_ranges = None
    self._query_sequence = None
    self._query_length = 0
    self._target_length = 0
    self._reference = None
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
  def get_reference(self):
    return self._reference

  def print_alignment(self,chunk_size=40):
    trantab = maketrans('01',' *')
    qseq = self.get_query_sequence()
    if self.get_strand() == '-': qseq = rc(qseq)
    prevt = self._alignment_ranges[0][0].start-1
    prevq = self._alignment_ranges[0][1].start-1
    qstrs = []
    tstrs = []
    qaln = ''
    taln = ''
    for r in self._alignment_ranges:
      t = r[0]
      q = r[1]
      diffq = q.start-prevq-1
      difft = t.start-prevt-1
      prevt = r[0].end
      prevq = r[1].end
      if difft > 68:
        qstrs.append(qaln)
        tstrs.append(taln)
        qaln = ''
        taln = ''
        difft = 0
        diffq = 0
      qaln += qseq[q.start-1:q.end+difft].upper()+'-'*diffq
      taln += self.get_reference()[t.chr][t.start-1:t.end+diffq].upper()+'-'*difft
    if len(qaln) > 0:
      qstrs.append(qaln)
      tstrs.append(taln)
    for i in range(len(qstrs)):
      print 'Exon '+str(i+1)
      mm = ''.join([str(int(qstrs[i][j]!=tstrs[i][j] and qstrs[i][j]!='-' and tstrs[i][j]!='-')) for j in range(len(qstrs[i]))]).translate(trantab)
      q = qstrs[i]
      t =  tstrs[i]
      for y in [[mm[x:x+chunk_size],q[x:x+chunk_size],t[x:x+chunk_size]] for x in range(0,len(mm),chunk_size)]:
        print '  '+y[0]
        print 'Q '+y[1]
        print 'T '+y[2]
        print ''
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
