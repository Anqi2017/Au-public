import sys
from Bio.Range import merge_ranges

# Classes to help stream biological data

# Works for any stream with a 
# 1. read_entry 
# function and also 
# 2. get_range
# function for each of the objects streamed
class LocusStream:
  def __init__(self,stream):
    self.stream = stream
    self.current_range = None
    firstobj = self.stream.read_entry()
    if not firstobj: return
    self.current_range = firstobj.get_range()
    self.current_range.set_payload([firstobj])

  def __iter__(self):
    return self
  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r
    
  def read_entry(self):
    if not self.current_range:
      return None
    output = None
    while True:
      e = self.stream.read_entry()
      if e:
        rng = e.get_range()
        if not rng: continue # continue if nonetype for range
        if rng.overlaps(self.current_range):
          self.current_range.get_payload().append(e)
          if self.current_range.end < rng.end: self.current_range.end = rng.end
        else: 
          output = self.current_range
          self.current_range = rng
          self.current_range.set_payload([e])
          break
      else:
        output = self.current_range
        self.current_range = None
        break
    return output

# Take an array streams
# Each element should be sorted by position
# Streams need to have this method:
# read_entry
# Each entry should have a get_range element
class MultiLocusStream:
  def __init__(self,streams):
    self.streams = streams
    self.buffers = []
    # seed the buffers
    for i in range(0,len(streams)):
      entry = self.streams[i].read_entry()
      self.buffers.append(entry)
    #self.set_current_range()

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r

  def read_entry(self):
    # Find our current lowest range
    output = []
    for i in self.buffers: output.append([])
    rngs = [x.get_range() for x in self.buffers if x]
    #print rngs
    if len(rngs) == 0: return None
    srngs = sorted(rngs,key=lambda x: (x.chr,x.start,x.end))
    mrngs = merge_ranges(srngs)
    current_range = mrngs[0]
    #print current_range.get_range_string()
    got_overlap = True
    while got_overlap == True:
      got_overlap = False
      for i in range(0,len(self.buffers)):
        if not self.buffers[i]: continue #end of this one
        v = self.buffers[i].get_range().cmp(current_range)
        if v==0:
          got_overlap = True
          if self.buffers[i].get_range().overlaps(current_range):
            current_range = current_range.merge(self.buffers[i].get_range())
          output[i].append(self.buffers[i])
          self.buffers[i] = self.streams[i].read_entry()
          #print str([len(x) for x in output])+"\t"+current_range.get_range_string()
        #print str(i)+":"+str(v)
    current_range.set_payload(output)
    return current_range

# Take an array streams
# Each element should be sorted by position
# Streams need to have this method:
# read_entry
# Each entry should have a get_range element
class MultiLocusStream2:
  def __init__(self,streams):
    self.streams = streams
    # self.buffers holds a list for each stream
    self.buffers = []
    self.current_range = None
    # seed the buffers
    for i in range(0,len(streams)):
      entry = self.streams[i].read_entry()
      self.buffers.append([entry])
    self.set_current_range()

  def set_current_range(self):
    lowest = self.get_lowest()
    if lowest == None: return # nothing to change on range since everything is looking done
    self.current_range = self.buffers[lowest][-1].get_range().copy()
    for i in range(0,len(self.buffers)):
      if i==lowest: continue
      if not self.buffers[i][-1]: continue
      if self.buffers[i][-1].get_range().overlaps(self.buffers[lowest][-1].get_range()):
        self.current_range = self.current_range.merge(self.buffers[i][-1].get_range())

  #Post: return the index of the lowest current buffer
  def get_lowest(self):
    #One type defined by type= argument
    #print [len(x) for x in self.buffers]
    #if self.current_range:
    #  print self.current_range.length()
    #  print self.current_range.get_range_string()
    nonzero = [i for i in range(0,len(self.buffers)) if self.buffers[i][-1]]
    vs = sorted(nonzero, key = lambda x: (self.buffers[x][-1].get_range().chr, self.buffers[x][-1].get_range().start, self.buffers[x][-1].get_range().end))
    if len(vs) == 0: return None
    return vs[0]

  #compare the current range to the end of everything
  def check_status(self):
    buffer_cache = [self.current_range.cmp(x[-1].get_range()) for x in self.buffers if x[-1]]
    if len(buffer_cache) ==0:
      return -1
    return max(buffer_cache)
  
  def fill_lowest(self):
    added_cache = set()
    while self.check_status() != -1:
      lowest = self.get_lowest()
      if lowest == None: 
        #print 'lowest '+str(lowest)
        return None# we are done here
      entry = self.streams[lowest].read_entry()
      if not entry:
        #perhaps end of data
        self.buffers[lowest].append(None)
        continue
        #return None
      rng = entry.get_range()
      self.buffers[lowest].append(entry)
      if rng.get_range_string() not in added_cache:
        if entry.get_range().overlaps(self.current_range):
          self.current_range = self.current_range.merge(entry.get_range())
        added_cache.add(rng.get_range_string())
    return 1

  # Return the current entry
  def read_entry(self):
    # fill the lowest
    status = self.fill_lowest()
    r = [[y for y in self.buffers[x][:-1]] for x in range(0,len(self.buffers))]
    # r is the return value
    # now reset the buffer
    self.buffers = [[x[-1]] for x in self.buffers]
    rng = self.current_range
    rng.set_payload(r)
    self.set_current_range() # new current range to work in
    #finished = True
    #for s in self.streams:
    #  e = s.read_entry()
    #  if e: finished = False
    #  v.append(e)
    #if finished: return None
    if max([len(x) for x in r]) == 0: return None
    #if not rng.get_payload(): return None
    #if len(rng.get_payload())==0: return None
    return rng

  def __iter__(self):
    return self

  def next(self):
    r = self.read_entry()
    if not r: raise StopIteration
    else:
      return r
    
