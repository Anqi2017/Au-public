import sys
from SequenceBasics import rc

class PairwiseAlignment:
  def __init__(self):
    self.s1 = None
    self.s2 = None
    self.a1 = None
    self.a2 = None
    self.score = None
  def load_sequences(self,s1,s2):
    self.s1 = s1
    self.s2 = s2
  def do_needleman_wunsch(self):
    [a1,a2,score] = needleman_wunsch(self.s1,self.s2)
    self.a1 = a1
    self.a2 = a2
    self.score = score
  def load_alignment(self,a1,a2):
    self.a1 = a1
    self.a2 = a2
  def count_errors(self):
    if not self.a1 or not self.a2:
      sys.stderr.write("You need to load or generate alignments\n")
      sys.exit()
    if len(self.a1) != len(self.a2):
      sys.stderr.write("Your alignments need to be equal length\n")
      sys.exit()
    insertions_into_s1 = 0
    deletions_from_s1 = 0
    mismatches = 0
    matches = 0
    for i in range(0,len(self.a1)):
      if self.a1[i] == '-':
        insertions_into_s1 += 1
      elif self.a2[i] == '-':
        deletions_from_s1 += 1
      elif self.a1[i] != self.a2[i]:
        mismatches += 1
      else:
        matches += 1
    return [matches, mismatches, insertions_into_s1, deletions_from_s1]

def Snw(c1,c2):
  if c1 == c2: return 10
  else: return -5
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

# The Alignment result from SmithWatermanAligner
class SmithWatermanAlignment:
  def __init__(self):
    return
  def set_alignment(self,gapopen,gapextend,match,mismatch,bidirectional,score,a1,a2,start_1,start_2,strand_1,strand_2,s1,s2):
    self.parameters = {}
    self.parameters['gapopen'] = gapopen
    self.parameters['gapextend'] = gapextend
    self.parameters['match'] = match
    self.parameters['mismatch'] = mismatch
    self.parameters['bidirectional'] = bidirectional
    self.score = score
    self.alignment_1 = a1
    self.alignment_2 = a2
    self.strand_1 = strand_1
    self.strand_2 = strand_2
    self.start_1 = start_1
    self.start_2 = start_2
    self.sequence_1 = s1
    self.sequence_2 = s2
    return
  def print_alignment(self):
    print str(self.parameters)
    print 'Score: ' +str(self.score)
    print  self.strand_1+" "+str(self.start_1)\
          +' '*max(0,(len(str(self.start_2))-len(str(self.start_1))))+' '\
          +self.alignment_1
    print  self.strand_2+" "+str(self.start_2)\
          +' '*max(0,(len(str(self.start_1))-len(str(self.start_2))))+' '\
          +self.alignment_2

class SmithWatermanAligner:
  def __init__(self):
    # Run parameters
    self.gapextend = -4
    self.gapopen = -10
    self.match = 10
    self.mismatch = -15
    self.bidirectional = True # Try both directions of s2 if true
    # User input sequences
    self.input_s1 = None
    self.input_s2 = None
    # Running variables
    self.M = None
    self.s1 = None
    self.s2 = None
    return

  def align(self):
    self.s1 = self.input_s1
    self.s2 = self.input_s2
    self.M = None
    outs1 = self.execute_sw_alignment2()
    outs1.append('+')
    outs1.append('+')
    if self.bidirectional:
      self.s1 = self.input_s1
      self.s2 = rc(self.input_s2)
      self.M = None
      outs2 = self.execute_sw_alignment2()
      outs2.append('+')
      outs2.append('-')
      if outs2[0] > outs1[0]:
        outs1 = outs2
    result = SmithWatermanAlignment()
    result.set_alignment(self.gapopen,self.gapextend,self.match,\
                         self.mismatch,self.bidirectional,outs1[0],outs1[1],\
                         outs1[2],outs1[3],outs1[4],outs1[5],outs1[6],\
                         self.input_s1,self.input_s2)
    return result

  def set_unidirectional(self):
    self.bidirectional = False
  def set_bidirectional(self):
    self.bidirectional = True

  def set_sequences(self,s1,s2):
    self.input_s1 = s1
    self.input_s2 = s2
    return


  # Pre: A 2D array
  # Post: Prints matrix to STDOUT
  # Modifies: STDOUT
  def print_matrix(self):
    M = self.M
    if len(M) == 0: return
    for m in range(0,len(M)):
      oline = ''
      for n in range(0,len(M[0])):
        oline = oline + ' ' + str(M[m][n])
      print oline

  # Pre: Two sequences and a 2D array
  def print_alignment_matrix(self):
    s1= self.s1
    s2 = self.s2
    M = self.M
    if len(M) == 0: return
    if len(M) != len(s2) or len(M[0]) != len(s1): return
    line1 = ' '
    for c in s1:
      line1 = line1 +  ' ' + str(c)
    print line1
    oline = ''
    for m in range(0,len(M)):
      oline = s2[m] + ''
      for n in range(0,len(M[0])):
        oline = oline + ' ' + str(M[m][n])
      print oline

  #Fetch the diagnal M value from the current coordinate
  def diag_score(self,i,j):
    if i-1 < 0 or j-1 < 0: return 0
    return self.M[i-1][j-1]

  # Return the score given the current characters
  def match_score(self,c1,c2):
    if c1 == c2: return self.match
    return self.mismatch

  # Pre matrix M and current position i, j
  # Post: An array of scores
  def row_scores(self,i,j):
    oscores = []
    if i == 0:
      oscores.append(0)
      return oscores
    bottom = 0
    if i-self.maxgap > 0: bottom = i - self.maxgap
    for m in range(bottom,i):
      k = i-m #distance
      oscores.append(self.M[i-k][j]+self.gapopen+(k-1)*self.gapextend)
    return oscores

  # Pre: Score matrix M and current position i j
  # Post: An array of scores
  def col_scores(self,i,j):
    oscores = []
    if j == 0:
      oscores.append(0)
      return oscores
    bottom = 0
    if j - self.maxgap > 0: bottom = j - self.maxgap
    for m in range(bottom,j):
      l = j - m #distance
      oscores.append(self.M[i][j-l]+self.gapopen+(l-1)*self.gapextend)
    oscores.append(self.M[i][j-1]+self.gapopen)
    return oscores

  #Make the M scoring matrix for the alignment
  # Pre: Two sequences
  # Post: M a matrix with the scores computed
  # Modifies: Sets/alters self.M
  def score_matrix(self):
    s1 = self.s1
    s2 = self.s2
    self.M = [[0 for x in range(0,len(s1))] for x in range(0,len(s2))] #initialize alignment matrix
    for m in range(0,len(self.M)):
      for n in range(0,len(self.M[0])):
        #print self.col_scores(m,n)
        #sys.exit()
        updated_value = max(self.diag_score(m,n) + self.match_score(s1[n],s2[m]),max(self.row_scores(m,n)),max(self.col_scores(m,n)),0)
        self.M[m][n] = updated_value

  # Pre: Matrix M
  # Post: list [i,j] with the best coordinate
  def matrix_max(self):
    M = self.M
    best = [0,0]
    bestval = 0
    for m in range(0,len(M)):
      for n in range(0,len(M[0])):
        if M[m][n] > bestval:
          best = [m,n]
          bestval = M[m][n]
    return best

  # Pre M scoring matrix, current coordinate
  # Post: next score and coordinate
  def next_coord(self,i,j):
    M = self.M
    rowval = 0
    if i -1 >= 0:
      rowval =M[i-1][j]
    colval = 0
    if j-1 >= 0:
      colval = M[i][j-1]
    diagval = 0
    if i-1 >= 0 and j-1 >= 0:
      diagval = M[i-1][j-1] 
    if diagval >= rowval and diagval >= colval:
      return [diagval,i-1,j-1]
    if rowval >= colval:
      return [rowval, i-1, j]
    return [colval,i,j-1]

  def execute_sw_alignment(self):
    s1 = self.s1
    s2 = self.s2
    self.score_matrix()
    M = self.M
    [i,j] = self.matrix_max()
    currentscore = M[i][j]
    maxscore = currentscore
    s1o = []
    s2o = []
    a1 = ''
    a2 = ''
    [isave, jsave] = [0,0]
    while currentscore > 0 and i >= 0 and j >= 0:
      [isave, jsave] = [i,j]
      [currentscore,inext,jnext] = self.next_coord(i,j)
      if inext == i: #skip one on s2
        s1o.insert(0,s1[j])
        s2o.insert(0,'-')
      elif jnext == j: #skip one on s1
        s1o.insert(0,'-')
        s2o.insert(0,s2[i])
      else:
        s1o.insert(0,s1[j])
        s2o.insert(0,s2[i])
      [i,j] = [inext, jnext]
    s1start = jsave+1
    s2start = isave+1
    a1 = ''.join(s1o)
    a2 = ''.join(s2o)
    return [maxscore, a1, a2, s1start, s2start]

  def execute_sw_alignment2(self):
    self.M = []
    #Initialize Matrix
    for i in range(0,len(self.s2)+1):
      self.M.append([])
      for j in range(0,len(self.s1)+1):
        self.M[i].append({})
        self.M[i][j]['score'] = 0
        self.M[i][j]['pointer'] = 'none'
    #Fill matrix
    max_i = 0
    max_j = 0
    max_score = 0
    for i in range(1,len(self.s2)+1):
      for j in range(1,len(self.s1)+1):
        diag_score = 0
        left_score = 0
        up_score = 0
        chr1 = self.s1[j-1]
        chr2 = self.s2[i-1]
        if chr1 == chr2:
          diag_score = self.M[i-1][j-1]['score']+self.match
        else:
          diag_score = self.M[i-1][j-1]['score']+self.mismatch
        if self.M[i-1][j]['pointer'] == 'up':
          upscore = self.M[i-1][j]['score']+self.gapextend
        else:   
          up_score = self.M[i-1][j]['score']+self.gapopen
        if self.M[i-1][j]['pointer'] == 'left':
          left_score = self.M[i][j-1]['score']+self.gapextend
        else:
          left_score = self.M[i][j-1]['score']+self.gapopen
        if diag_score <= 0 and up_score <= 0 and left_score <= 0:
          self.M[i][j]['score'] = 0
          self.M[i][j]['pointer'] = 'none'
          continue
        if diag_score >= up_score and diag_score >= left_score:
          self.M[i][j]['score'] = diag_score
          self.M[i][j]['pointer'] = "diagonal"
        elif up_score >= left_score:
          self.M[i][j]['score'] = up_score
          self.M[i][j]['pointer'] = "up"
        else:
          self.M[i][j]['score'] = left_score
          self.M[i][j]['pointer'] = 'left'
        if self.M[i][j]['score'] > max_score:
          max_i = i
          max_j = j
          max_score = self.M[i][j]['score']
    a1 = ''
    a2 = ''
    j = max_j
    i = max_i
    while True:
      if self.M[i][j]['pointer'] == 'none': break
      if self.M[i][j]['pointer'] == 'diagonal':
        a1 = self.s1[j-1]+a1
        a2 = self.s2[i-1]+a2
        i-=1
        j-=1
      elif self.M[i][j]['pointer'] == 'left':
        a1 = self.s1[j-1]+a1
        a2 = '-'+ a2 
        j-=1
      elif self.M[i][j]['pointer'] == 'up':
        a1 = '-'+a1
        a2 = self.s2[i-1]+a2
        i-=1
    return [max_score, a1, a2, j, i]
