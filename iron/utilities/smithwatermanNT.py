#!/usr/bin/python
import sys

### smithwaterman.py ###
# A very basic Smith-Waterman alignment
# Input: two sequences
# Output: Print a local alignment with the starting nucleotide of each sequence before each sequence
# Modifies: STDOUT

maxgap = 10
gapopen = -5
gapextend = -5 
match = 10
mismatch = -5

#### print_matrix ###
# Input: A 2D array
# Output:  Print matrix to STDOUT
# Modifies:  STDOUT
def print_matrix(M):
  if(len(M)==0):
    return
  mlen = len(M)
  nlen = len(M[0])
  for m in range(0,mlen):
    oline = ''
    for n in range(0,nlen):
      oline = oline + ' ' + str(M[m][n])
    print oline

#### print_alignment_matrix ###
# Input: Two sequences, and a 2D array
# Output:  Print alignment matrix to STDOUT
# Modifies:  STDOUT
def print_alignment_matrix(M,s1,s2):
  if(len(M)==0): return
  if(len(M)!=len(s2) or len(M[0]) != len(s1)): return
  mlen = len(M)
  nlen = len(M[0])
  line1 = ' '
  for c in s1:
    line1 = line1 + ' ' + str(c)
  print line1
  for m in range(0,mlen):
    oline = s2[m] + ''
    for n in range(0,nlen):
      oline = oline + ' ' + str(M[m][n])
    print oline

#### diag_score ####
# Fetch the diagnal H value from the current coordinate
# Input:  Matrix H and current coordiante row i and column j
# Output: The Hi-1,j-1 value
# Modifies: None
def diag_score(H,i,j):
  if(i-1 < 0 or j-1 < 0): return 0
  return H[i-1][j-1]

#### match_score ####
# Return the score given the current characters
# Input:  Two characters c1 c2
# Output: The score
# Modifies: None
def match_score(c1,c2):
  if(c1 == c2): return match
  return mismatch

#### row_scores ####
# Return the scores for the gap going up the row
# Input:  Score matrix H and current position i j
# Output: an array of scores
# Modifies: None
def row_scores(H,i,j):
  oscores = list()
  if(i==0):
    oscores.append(0)
    return oscores
  bottom = 0
  if i-maxgap > 0: bottom = i-maxgap
  for m in range(bottom,i):
    k=i-m #distance
    oscores.append(H[i-k][j]+gapopen+(k-1)*gapextend)
  return oscores

#### col_scores ####
# Return the scores for the gap going across the columnb
# Input:  Score matrix H and current position i j
# Output: an array of scores
# Modifies: None
def col_scores(H,i,j):
  oscores = list()
  if(j==0):
    oscores.append(0)
    return oscores
  bottom = 0
  if j-maxgap > 0: bottom = j-maxgap
  for m in range(bottom,j):
    l=j-m #distance
    oscores.append(H[i][j-l]+gapopen+(l-1)*gapextend)
  return oscores

#### score_matrix ###
# Make the H scoring matrix for the alginment
# Input: Two sequences
# Output:  H a matrix with with the scores computed
# Modifies:  STDOUT
def score_matrix(s1,s2):
  H = [[0 for x in range(0,len(s1))] for x in range(0,len(s2))] #initialize alignment matrix
  mlen = len(H)
  nlen = len(H[0])
  for m in range(0,mlen):
    for n in range(0,nlen):
      H[m][n] = max(diag_score(H,m,n) + match_score(s1[n],s2[m]),max(row_scores(H,m,n)),max(col_scores(H,m,n)),0)
      #print_alignment_matrix(H,s1,s2)
      #print ''
  return H

#### matrix_max #####
# return the coordinate of the max value
# Input: takes a matrix H
# Output: list [i,j] with the best coordiante
def matrix_max(H):
  mlen = len(H)
  nlen = len(H[0])
  best = [0,0]
  bestval = 0
  for m in range(0,mlen):
    for n in range(0,nlen):
      if H[m][n] > bestval: 
        best = [m , n]
        bestval = H[m][n]
  return best

#### next_coord ###
# Print the next coordinate to go to and its score
# Input: H scoring matrix, and current coordinate i j
# Output:  the next score and coordiante inext jnext
# Modifies:  None
def next_coord(H,i,j):
  rowval = 0
  if(i-1 >= 0):
    rowval = H[i-1][j]
  colval = 0
  if(j-1 >= 0):
    colval = H[i][j-1]
  diagval = 0
  if(i-1 >=0 and j-1 >= 0):
    diagval = H[i-1][j-1]
  if(diagval >= rowval and diagval >= colval):
    return [diagval,i-1,j-1]
  if(rowval >= colval):
    return [rowval,i-1,j]
  return [colval,i,j-1]


#### get_local_alignment ###
# Print the local alignment given the scoring matrix
# Input: H scoring matrix, sequences s1 and s2
# Output:  A best local alignment between the two sequences, returns the max alignment score, and two strings that are the alignment lines, and the start indecies for the two sequences, and a descriptor of how s2 differs from s1
# Modifies:  none
def get_local_alignment(H,s1,s2):
  mlen = len(H)
  nlen = len(H[0])
  [i,j] = matrix_max(H)
  currentscore = H[i][j]
  maxscore = currentscore
  s1o = list()
  s2o = list()
  a1 = ''
  a2 = ''
  [isave,jsave] = [0,0]
  while(currentscore > 0 and i >= 0 and j >= 0):
    [isave, jsave] = [i,j]
    [currentscore,inext,jnext] = next_coord(H,i,j)
    if(inext==i): # skip one on s2
      s1o.insert(0,s1[j])
      s2o.insert(0,'-')
    elif(jnext==j): #skip one on s1
      s1o.insert(0,'-')
      s2o.insert(0,s2[i])
    else:
      s1o.insert(0,s1[j])
      s2o.insert(0,s2[i])
    [i,j] = [inext,jnext]
  s1start = jsave+1
  s2start = isave+1
  a1 = ''.join(s1o)
  a2 = ''.join(s2o)
  return [maxscore, a1, a2,s1start,s2start]

s1 = sys.argv[1]
s2 = sys.argv[2]

M = [[0 for x in range(0,len(s1))] for x in range(0,len(s2))] #initialize alignment matrix

H = score_matrix(s1,s2)

#print_alignment_matrix(H,s1,s2)
[maxscore,s1align,s2align,s1coord,s2coord] = get_local_alignment(H,s1,s2)
print str(maxscore) + "\t" + str(s1coord) + "\t" + s1align + "\t" + str(s2coord) + "\t" + s2align
