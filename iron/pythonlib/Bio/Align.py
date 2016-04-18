import re, sys
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
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
  def get_query_sequence(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
    return self._query_sequence
  def get_target_length(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
    return self._target_length
  def get_strand(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
    return self._query_direction
  def get_reference(self):
    sys.stderr.write("ERROR: needs overridden\n")
    sys.stderr.exit()
    return self._reference

  
  # Process the alignment to get information like
  # the alignment strings for each exon
  def _get_alignment_strings(self,min_intron_size=68):
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
      if difft >= min_intron_size:
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
    return [qstrs,tstrs]

  def _analyze_alignment(self,min_intron_size=68):
    [qstrs,tstrs] = self._get_alignment_strings(min_intron_size=68)
    matches = sum([x[0].length() for x in self._alignment_ranges]) 
    misMatches = 0
    for i in range(len(qstrs)):
      misMatches += sum([int(qstrs[i][j]!=tstrs[i][j] and qstrs[i][j]!='-' and tstrs[i][j]!='-' and tstrs[i][j]!='N') for j in range(len(qstrs[i]))])
    nCount = sum([len([y for y in list(x) if y == 'N'])  for x in tstrs])
    qNumInsert = sum([len(re.findall('[-]+',x)) for x in tstrs])
    qBaseInsert = sum([len(re.findall('[-]',x)) for x in tstrs])
    tNumInsert = sum([len(re.findall('[-]+',x)) for x in qstrs])
    tBaseInsert = sum([len(re.findall('[-]',x)) for x in qstrs])
    matches = matches - misMatches - nCount
    return {'matches':matches,\
            'misMatches':misMatches,\
            'nCount':nCount,\
            'qNumInsert':qNumInsert,\
            'qBaseInsert':qBaseInsert,\
            'tNumInsert':tNumInsert,\
            'tBaseInsert':tBaseInsert}

  def print_alignment(self,chunk_size=40,min_intron_size=68):
    trantab = maketrans('01',' *')
    [qstrs,tstrs] = self._get_alignment_strings(min_intron_size=68)
    for i in range(len(qstrs)):
      print 'Exon '+str(i+1)
      mm = ''.join([str(int(qstrs[i][j]!=tstrs[i][j] and qstrs[i][j]!='-' and tstrs[i][j]!='-' and tstrs[i][j]!='N')) for j in range(len(qstrs[i]))]).translate(trantab)
      q = qstrs[i]
      t =  tstrs[i]
      for y in [[mm[x:x+chunk_size],q[x:x+chunk_size],t[x:x+chunk_size]] for x in range(0,len(mm),chunk_size)]:
        print '  '+y[0]
        print 'Q '+y[1]
        print 'T '+y[2]
        print ''

  # These methods may be overrriden by an alignment type to just return themself
  # clearly this should be over written by the PSL type to just give itself
  def get_PSL(self,min_intron_size=68):
    from Bio.Format.PSL2 import PSL
    if not self._alignment_ranges: return None
    matches = sum([x[0].length() for x in self._alignment_ranges]) # 1. Matches - Number of matching bases that aren't repeats
    misMatches = 0 # 2. Mismatches - Number of baess that don't match
    repMatches = 0 # 3. repMatches - Number of matching baess that are part of repeats
    nCount = 0 # 4. nCount - Number of 'N' bases
    qNumInsert = 0 # 5. qNumInsert - Number of inserts in query
    qBaseInsert = 0 # 6. qBaseInsert - Number of bases inserted into query
    tNumInsert = 0 # 7. Number of inserts in target
    tBaseInsert = 0 # 8. Number of bases inserted into target
    sub = self.get_query_sequence()
    ref = self.get_reference()
    if sub and ref:
      v = self._analyze_alignment(min_intron_size=min_intron_size)
      matches = v['matches']
      misMatches = v['misMatches'] # 2. Mismatches - Number of baess that don't match
      nCount = v['nCount'] # 4. nCount - Number of 'N' bases
      qNumInsert = v['qNumInsert'] # 5. qNumInsert - Number of inserts in query
      qBaseInsert = v['qBaseInsert'] # 6. qBaseInsert - Number of bases inserted into query
      tNumInsert = v['tNumInsert'] # 7. Number of inserts in target
      tBaseInsert = v['tBaseInsert'] # 8. Number of bases inserted into target
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
    return PSL(psl_string,query_sequence=self.get_query_sequence(),reference=self.get_reference())

  #clearly this should be overwritten by the SAM class to give itself
  def get_SAM(self):
    return
