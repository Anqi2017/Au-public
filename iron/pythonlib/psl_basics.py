# psl_basics.py
#
# A module for holding functions for handling psl files
# 
# These include:
# read_psl_entry() - read a line from a psl file into a dictionary
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

import re

# pre: one psl entry
# post: one genepred line of the query's location on the target
#       we will use the genepred format with an aditional gene name field
def convert_entry_to_genepred_line(psl):
  vals = [psl['qName'], psl['qName'], psl['tName'], psl['strand'], str(psl['tStart']), str(psl['tEnd']), str(psl['tStart']), str(psl['tEnd']), str(psl['blockCount'])]
  starts = []
  ends = []
  for i in range(0,psl['blockCount']):
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

# pre: a line from a psl file
# post: a dictionary with terms assigned for each variable in the psl line
#       includes an additional 'qStarts_actual' keyed array that contains the
#         start locations for the query in the same coordinate system as the
#         it would have on the positive strand to simplify comparisions of 
#         coordinate overlaps between multiple hits of the same query regardless
#         of direction.
# modifies: none
def read_psl_entry(line):
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
