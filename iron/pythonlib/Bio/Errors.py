from Bio.Sequence import rc
import sys
### Error Analysis ####
# I am to describe errors at several levels
# 
# Errors in the query sequence
# 
# 1. Is a query base an error or not?
#    Probability - Sometimes it can be ambiguous which base is in error
# 
# 2. What is the basic type of error?
#    Mismatch
#    Insertion
#      Total insertion
#      Homopolymer insertion
#    Deletion
#      Total deletion
#        Before
#        After
#      Homopolymer deletion
#   sum of probabilities should add up to 1.)
#
# 3. What is the more specific error?
#    Mismatch type
#    insertion/deletion - Base, Length
#
valid_types = set(['match','mismatch','total_insertion','total_deletion','homopolymer_insertion','homopolymer_deletion'])


# Error to describe a single base in the query sequence
# A deletion could be either preceeding or following a base
class QueryBaseError():
  global valid_types
  def __init__(self):
    self._error_probability = 0
    self._total_deletion_before_probability = 0
    self._total_deletion_after_probability = 0
    self._total_deletion_before_details = {'base':'N','length':0}
    self._total_deletion_after_details = {'base':'N','length':0}
    self._deletion_before_probability = 0
    self._deletion_after_probability = 0
    self._type = None
    self._insert_details = None
    self._deletion_details = None
    self._single_base_details = None
    self._deletion_before_details = {'base':'N','length':0}
    self._deletion_after_details = {'base':'N','length':0}
  def get_adjusted_error_count(self):
    p1 = self._deletion_before_probability*self._deletion_before_details['length']+self._deletion_after_probability*self._deletion_after_details['length']+self._error_probability
    p2 = self._total_deletion_before_probability*self._total_deletion_before_details['length']+self._total_deletion_after_probability*self._total_deletion_after_details['length']
    return p1+p2
  def set_error_probability(self,prob):
    self._error_probability = float(prob)
  def get_error_probability(self):
    return self._error_probability
  def set_type(self,type):
    if type not in valid_types:
      sys.stderr.write("ERROR: type is not a valid type")
      sys.exit()
    self._type = type
  def set_total_deletion_before(self,p,base,dlen):
    self._total_deletion_before_details = {'base':base,'length':dlen}
    self._total_deletion_before_probability = p
  def set_total_deletion_after(self,p,base,dlen):
    self._total_deletion_after_details = {'base':base,'length':dlen}
    self._total_deletion_after_probability = p
  def get_type(self):
    return self._type
  def set_insert_details(self,base,ilen):
    self._insert_details = {'base':base,'length':ilen}
  def set_single_base_details(self,tbase,qbase):
    self._single_base_details = {'target':tbase,'query':qbase}
  def set_deletion_before(self,prob,base,dlen):
    self._deletion_before_probability = float(prob)
    self._deletion_details = {'base':base,'length':dlen}
  def set_deletion_after(self,prob,base,dlen):
    self._deletion_after_probability = float(prob)
    self._deletion_details = {'base':base,'length':dlen}
  def __str__(self):
    ostr = "Type: "+self._type+"\n"
    ostr += "P(error): "+str(self._error_probability)+"\n"
    if self.get_type() == 'total_insertion' or self.get_type() == 'homopolymer_insertion': 
      ostr += "Homopolymer Insertion length: "+str(self._insert_details['length'])+"\n"
    if self._deletion_before_probability > 0:
      ostr += "P(err-HP-del-before): "+str(self._deletion_before_probability)+"\n"
    if self._deletion_after_probability > 0:
      ostr += "P(err-HP-del-after): "+str(self._deletion_after_probability)+"\n"
    if self.get_type() == 'total_deletion' or self.get_type() == 'homopolymer_deletion':
      ostr += "Homopolymer Deletion length: "+str(self._deletion_details['length'])+"\n"
    if self._total_deletion_before_probability > 0:
      ostr += "P(err-total-del-before): "+str(self._total_deletion_before_probability)+"\n"
      ostr += "Total deletion before length: "+str(self._total_deletion_before_details['length'])+"\n"
    if self._total_deletion_after_probability > 0:
      ostr += "P(err-total-del-after): "+str(self._total_deletion_after_probability)+"\n"
      ostr += "Total deletion after length: "+str(self._total_deletion_after_details['length'])+"\n"
    return ostr

class AlignmentErrors:
  # Pre: Take an alignment between a target and query
  #      Uses get_strand from alignment to orient the query
  #      All results are on the positive strand of the query
  #      (meaning may be the reverse complement of target if negative)
  def __init__(self,alignment,min_intron_size=68):
    #self._alns = []
    self._min_intron_size=min_intron_size
    self._aligned_query = None
    self._hpas = []
    self._has_quality = False # can be changed when add_alignment uses one that has quality
    self._alignment = alignment
    self._quality_distro = None # gets set by analyze_quality
    self._deletion_type = None

    astrings = self._alignment.get_alignment_strings(min_intron_size=self._min_intron_size)
    if self._alignment.get_query_quality(): self._has_quality = True
    if len(astrings) == 0: return None
    alns = []
    for i in range(len(astrings[0])):
      if self._alignment.get_strand() == '+':
        alns.append({'query':astrings[0][i],'target':astrings[1][i],'query_quality':astrings[2][i]})
      else:
        alns.insert(0,{'query':rc(astrings[0][i]),'target':rc(astrings[1][i]),'query_quality':astrings[2][i][::-1]})
    #if self._alignment.get_strand() == '-':
    #  alns = alns[::-1]
    #get homopolymer alignments
    self._hpas = self._misalign_split(alns) # split alignment into homopolymer groups
    self._query_hpas = []
    self._target_hpas = []
    qi = 0
    for i in range(len(self._hpas)):
      prev = None
      if i > 0: prev = self._hpas[i-1]
      foll = None
      if i + 1 < len(self._hpas): foll = self._hpas[i+1]
      qlen = len(self._hpas[i].get_query())
      for j in range(0,qlen):
        self._query_hpas.append({'hpa':self._hpas[i],'pos':j,'prev-hpa':prev,'next-hpa':foll})
      qi+=qlen
    ti = 0
    for i in range(len(self._hpas)):
      prev = None
      if i > 0: prev = self._hpas[i-1]
      foll = None
      if i + 1 < len(self._hpas): foll = self._hpas[i+1]
      tlen = len(self._hpas[i].get_target())
      for j in range(0,tlen):
        self._target_hpas.append({'hpa':self._hpas[i],'pos':j,'prev-hpa':prev,'next-hpa':foll})
      ti+=tlen
    #print '----'
    #print self.get_query_sequence()
    #print self.get_target_sequence()
    #for i in range(0,len(self._query_hpas)):
    #  print self.get_query_error(i)
      
  def get_query_errors(self):
    v = []
    for i in range(len(self._query_hpas)):
      v.append(self.get_query_error(i))
    return v

  # Pre:  given an index in the aligned query
  # Post: return the error description for that base
  def get_query_error(self,i):
    x = self._query_hpas[i]
    h = x['hpa']
    pos = x['pos']
    prob = 0
    be = QueryBaseError()
    if len(h.get_query()) == len(h.get_target()) and \
      h.get_query()==h.get_target():
      be.set_type('match')
      be.set_single_base_details(h.get_query()[0],h.get_target()[0])
      # they perfectly match we won't need to set an error for the base
    elif h.get_nt() == '*': #mismatch
      be.set_type('mismatch')
      be.set_error_probability(1)
      be.set_single_base_details(h.get_query()[0],h.get_target()[0])
    elif len(h.get_target()) == 0: # total insertion
      be.set_type('total_insertion')
      be.set_error_probability(1)
      be.set_insert_details(h.get_query()[0],len(h.get_query()))
    elif len(h.get_query()) > len(h.get_target()): # homopolymer insertion ... might be the base that shouldn't have been inserted
      be.set_type('homopolymer_insertion')
      be.set_error_probability(float(len(h.get_target()))/float(len(h.get_query())))
      be.set_insert_details(h.get_query()[0],len(h.get_query())-len(h.get_target()))
    elif len(h.get_target()) > len(h.get_query()): # homopolymer deletion ... might be the base that shouldn't have been inserted
      be.set_type('match')
      be.set_single_base_details(h.get_query()[0],h.get_target()[0])
      be.set_type('homopolymer_deletion')
      dlen = len(h.get_target())-len(h.get_query())
      p = 1/float(len(h.get_query())+1)
      be.set_deletion_before(p,h.get_query()[0],dlen)
      be.set_deletion_after(p,h.get_query()[0],dlen)
      #be.set_error_probability(float(len(h.get_target()))/float(len(h.get_query())))
    #print pos
    if i != 0 and pos == 0: # check for a total deletion before
      prev = x['prev-hpa']
      if len(prev.get_query()) == 0: # total deletion
        dlen = len(prev.get_target())
        be.set_total_deletion_before(0.5,prev.get_target()[0],dlen)
    if i != len(self._query_hpas)-1 and pos == len(h.get_query())-1: # check for a total deletion before
      if x['next-hpa']:
        foll = x['next-hpa']
        if len(foll.get_query()) == 0: # total deletion
          dlen = len(foll.get_target())
          be.set_total_deletion_after(0.5,foll.get_target()[0],dlen)
    #else: print h
    return be
    #elif len(h.target()) > len(h.query()): # homopolymer deletion ... might be the base that shouldn't have been inserted
    #  be.set_error_probability(float(len(h.target()))/float(len(h.query())))

  def get_query_sequence(self):
    return ''.join([x['hpa'].get_query()[0] for x in self._query_hpas])
  def get_target_sequence(self):
    return ''.join([x['hpa'].get_target()[0] for x in self._target_hpas])

  # Go through HPAGroups and store the distro of ordinal values of quality scores
  def analyze_quality(self):
    res = {}
    for h in self._hpas:
      if h.type() not in res: res[h.type()]={}
      for c in h.get_quality():
        if c not in res[h.type()]: res[h.type()][c] = 0
        res[h.type()][c]+=1
    self._quality_distro = res
  def get_quality_report_string(self):
    if not self._quality_distro:
      self.analyze_quality()
    ostr = ""
    for type in sorted(self._quality_distro.keys()):
      total = sum([ord(x)*self._quality_distro[type][x] for x in self._quality_distro[type]])
      cnt = sum([self._quality_distro[type][x] for x in self._quality_distro[type]])
      if cnt == 0: continue
      print 'type: '+type+' '+str(cnt)+' '+str(float(total)/float(cnt))
    return ostr

  def has_quality(self):
    return self._has_quality
  #Pre: alignment strings have been set so for each exon we have
  #     query, target and query_quality
  #     _has_quality will specify whether or not the quality is meaningful
  def _misalign_split(self,alns):
    total = []
    z = 0
    for x in alns:
      z += 1
      exon_num = z
      if self._alignment.get_strand() == '-':
        exon_num = (len(alns)-z)+1
      buffer = {'query':x['query'][0],'target':x['target'][0],'query_quality':x['query_quality'][0],'exon':exon_num}
      if buffer['query'] == '-': buffer['nt'] = buffer['target']
      elif buffer['target'] == '-': buffer['nt'] = buffer['query']
      elif buffer['query'] == buffer['target']: buffer['nt'] = buffer['query']
      elif buffer['query'] != buffer['target']: buffer['nt'] = '*'
      else:
	sys.stderr.write("WARNING unkonwn case\n")
      for i in range(1,len(x['query'])):
        qchar = x['query'][i]
        tchar = x['target'][i]
        qualchar = x['query_quality'][i]
        if qchar != tchar and (qchar != '-' and tchar != '-'):
          #classic mismatch
          #print 'mismatch'
          #print buffer
          total.append(buffer)
          buffer = {'query':qchar,'target':tchar,'query_quality':qualchar,'exon':exon_num}
          buffer['nt'] = '*'
        elif qchar == buffer['nt'] or tchar == buffer['nt']:
          # its a homopolymer match
          buffer['query'] += qchar

          buffer['target'] += tchar
          buffer['query_quality'] += qualchar
          #print 'homopoly'
        else:
          #print 'new thing'
          #print buffer
          total.append(buffer)
          buffer = {'query':qchar,'target':tchar,'query_quality':qualchar,'exon':exon_num}
          if qchar == '-': buffer['nt'] = tchar
          else:  buffer['nt'] = qchar
      total.append(buffer)
    result = [AlignmentErrors.HPAGroup(self,y) for y in total]
    return result

  #Homopolymer alignment group
  class HPAGroup:
    # takes a chunk of homopolymer alignment
    # as a dictionary with 'query' and 'target' sequences set
    # query should always be positive strand
    def __init__(self,parent,mydict):
      self._error_profile = parent
      self._data = mydict
      self._qseq = self._data['query'].replace('-','')
      self._tseq = self._data['target'].replace('-','')
      self._nt = self._data['nt'] # the nulceotide or * for mismatch
      self._qquality = self._data['query_quality'].replace('\0','')
      self._exon_number = self._data['exon']
      self._type = None
      if self._qseq == self._tseq:
        self._type = 'match'
      ### handle mismatches
      elif self._nt == '*':
        self._type = 'mismatch'
        self._code = self._tseq+'>'+self._qseq
      # Total deletion
      elif len(self._qseq) == 0:
        self._type = 'total_deletion'
      # Total insert
      elif len(self._tseq) == 0:
        self._type = 'total_insertion'
      elif len(self._qseq) < len(self._tseq):
        self._type = 'homopolymer_deletion'
      elif len(self._qseq) > len(self._tseq):
        self._type = 'homopolymer_insertion'
      else:
	sys.stderr.write("ERROR unsupported type\n")
        sys.exit()
    def get_nt(self): return self._data['nt']
    # always + strand
    def get_query(self):  return self._qseq
    # could be + or - strand
    def get_target(self):  return self._tseq
    def get_exon(self):  return self._exon_number
    def get_length(self):
      return {'query':len(self._qseq),'target':len(sef._tseq)}
    def __str__(self):
      return self.get_string()
    def get_string(self):
      ostr = ''
      ostr += 'Target:  '+self._tseq+"\n"
      ostr += 'Query:   '+self._qseq+"\n"
      if self._error_profile.has_quality(): ostr += 'Quality: '+self._qquality+"\n"
      ostr += 'Type: '+str(self._type)+"\n"
      return ostr
    def has_quality(self):
      return self._error_profile.has_quality()
    def get_quality(self):
      if not self.has_quality(): return False
      return self._qquality
    def type(self):
      return self._type


