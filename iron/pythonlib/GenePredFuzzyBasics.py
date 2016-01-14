# Class for finding general genepred matches
import sys
from GenePredBasics import GenePredEntry
from RangeBasics import Bed

# Similar to genepred but no requirement for exact borders
class FuzzyGenePred:
  #set use_dir true if you want to use direction and make it direction specific
  #set full_length true if you want to only bring together long reads that are full length matches
  def __init__(self,ingpd=None,use_dir=False,jun_tol=0,full_length=False,\
               do_add_single_exon=True,single_exon_minimum_length=200,\
               single_exon_minimum_overlap_fraction=0.8,\
               single_exon_minimum_overlap_bases=1,\
               single_exon_maximum_endpoint_distance=1000): 
    self.fuzzy_junctions = []
    self.gpds = [] #contributing member genepreds  
    self.dir = None
    self.use_dir = use_dir
    self.start = None
    self.end = None
    #Important parameter for how much tolerance at the junction site.
    self.junction_tolerance = jun_tol
    #Not fully implemented.  Do we require a full length match
    self.full_length = False
    # Define thresholds for overlapping single exons
    self.do_add_single_exon = do_add_single_exon
    self.single_exon_minimum_length = single_exon_minimum_length
    self.single_exon_minimum_overlap_fraction = single_exon_minimum_overlap_fraction #reciprocal ... must be this fraction or more on both
    self.single_exon_minimum_overlap_bases = single_exon_minimum_overlap_bases #minimum number of bases
    self.single_exon_maximum_endpoint_distance = single_exon_maximum_endpoint_distance
    if ingpd:
      self.add_gpd(ingpd)

  def exon_count(self):
    return len(self.fuzzy_junctions)+1
  def gpd_count(self):
    return len(self.gpds)
  #This is an inspection tool for a fuzzy gpd
  def get_info_string(self):
    ostr = ''
    ostr += "== FUZZY GENEPRED INFO =="+"\n"
    ostr += str(len(self.gpds))+' total GPDs'+"\n"
    totalbounds = Bed(self.start.chr,self.start.start,self.end.end,self.start.direction)
    ostr += totalbounds.get_range_string()+" total bounds\n";
    ostr += '---- start ----'+"\n"
    ostr += str(len(self.start.get_payload()))+ " reads supporting start"+"\n"
    ostr += '  '+str(mean(self.start.get_payload()))+' mean'+"\n"
    ostr += '  '+str(mode(self.start.get_payload()))+' mode'+"\n"
    ostr += '  '+self.start.get_range_string()+" start range\n"
    ostr += '---- end ----'+"\n"
    ostr += str(len(self.end.get_payload()))+ " reads supporting end"+"\n"
    ostr += '  '+str(mean(self.end.get_payload()))+' mean'+"\n"
    ostr += '  '+str(mode(self.end.get_payload()))+' mode'+"\n"
    ostr += '  '+self.end.get_range_string()+" end range\n"
    ostr += '---- junctions ----'+"\n"
    ostr += str(len(self.fuzzy_junctions))+' total fuzzy junctions'+"\n"
    cnt = 0
    for j in self.fuzzy_junctions:
      cnt += 1
      ostr += '  '+str(cnt)+'. '+str(mode(j.left.get_payload()['junc']))+" ^ "+str(mode(j.right.get_payload()['junc']))+"\n"
      ostr += "     "+j.left.get_range_string()+" ^ "+j.right.get_range_string()+"\n"
      ostr += "     "+str(len(j.left.get_payload()['junc']))+" read support" + "\n"
      if j.left.get_payload()['start']:
        ostr += "       "+"---starts----"+"\n"
        ostr += "       "+str(len(j.left.get_payload()['start'].get_payload()))+" starts at "+j.left.get_payload()['start'].get_range_string()+"\n"
      if j.right.get_payload()['end']:
        ostr += "       "+"---ends----"+"\n"
        ostr += "       "+str(len(j.right.get_payload()['end'].get_payload()))+" ends at "+j.right.get_payload()['end'].get_range_string()+"\n"
    return ostr

  #Add a new gpd return true if successful
  #Return false if it didn't work
  def add_gpd(self,ingpd):
    chr = ingpd.value('chrom')
    if len(self.gpds)==0:  # first one
      self.read_first(ingpd)
      return True
    # more difficult situation where we must try to combine
    # See if it can match first before actually adding stuff to it
    #if self. 
    newfuz = FuzzyGenePred(ingpd,use_dir=self.use_dir,full_length=self.full_length,jun_tol=self.junction_tolerance)
    self.add_fuzzy_gpd(newfuz)

  def add_fuzzy_gpd(self,fuz2):
    # see if we can add this fuzzy gpd to another
    # We treat single exon genes seprately so if only one of them is
    # single exon we can't compare them
    if len(fuz2.fuzzy_junctions) == 0 and len(self.fuzzy_junctions) != 0:
      return False
    if len(fuz2.fuzzy_junctions) != 0 and len(self.fuzzy_junctions) == 0:
      return False
    # Lets work combine the single exon step and exit
    if len(fuz2.fuzzy_junctions) == 0 and len(self.fuzzy_junctions) == 0:
      return self.do_add_single_exon_fuzzy_gpd(fuz2)
    #1. First we need perfect junctions for a run of them
    if not self.compatible_overlap(fuz2): return False
    # If they are both single exon genes we can just put them together
    if len(self.fuzzy_junctions)==0 and len(fuz2.fuzzy_junctions)==0:
      for s in fuz2.start.get_payload():
        self.start.get_payload().append(s)
      for e in fuz2.end.get_payload():
        self.end.get_payload().append(e)
      return True
    # For now don't add them if one is single exon
    if len(self.fuzzy_junctions)==0 or len(fuz2.fuzzy_junctions)==0:
      return False
    # If we are still here we know we can add the two of them together
    #print mode(fuz2.fuzzy_junctions[0].left.get_payload()['junc'])
    #print mode(self.fuzzy_junctions[0].left.get_payload()['junc'])
    # If they have the same starting junction we can add their starting points together
    if self.fuzzy_junctions[0].overlaps(fuz2.fuzzy_junctions[0],self.junction_tolerance):
      #print 'samestart'    
      for s in fuz2.start.get_payload():
        self.start.get_payload().append(s)
    # Check if the other one is new start
    elif mode(fuz2.fuzzy_junctions[0].left.get_payload()['junc']) < mode(self.fuzzy_junctions[0].left.get_payload()['junc']):
      #print "2 start"
      self.start = fuz2.start
    elif mode(fuz2.fuzzy_junctions[0].left.get_payload()['junc']) > mode(self.fuzzy_junctions[0].left.get_payload()['junc']):
      True
    #  #print "1 start"
    #  #we're good to go
    else:
      sys.stderr.write("WARNING: strange start case abort merge\n")
      return False
    # lets work the ends now
    if self.fuzzy_junctions[-1].overlaps(fuz2.fuzzy_junctions[-1],self.junction_tolerance):
      #print 'sameend'    
      for e in fuz2.end.get_payload():
        self.end.get_payload().append(e)
    # Check if the other one is new start
    elif mode(fuz2.fuzzy_junctions[-1].right.get_payload()['junc']) > mode(self.fuzzy_junctions[-1].right.get_payload()['junc']):
      #print "2 end"
      self.end = fuz2.end
    elif mode(fuz2.fuzzy_junctions[-1].right.get_payload()['junc']) < mode(self.fuzzy_junctions[-1].right.get_payload()['junc']):
      True
    #  #print "1 end"
    #  #we're good to go
    else:
      sys.stderr.write("WARNING: strange end case abort merge\n")
      u1= mode(self.fuzzy_junctions[-1].left.get_payload()['junc'])
      u2= mode(fuz2.fuzzy_junctions[-1].left.get_payload()['junc'])
      v1= mode(self.fuzzy_junctions[-1].right.get_payload()['junc'])
      v2= mode(fuz2.fuzzy_junctions[-1].right.get_payload()['junc'])
      sys.stderr.write(str(u1)+"\t"+str(u2)+"\n")
      sys.stderr.write(str(v1)+"\t"+str(v2)+"\n")
      return False
    # now the starts and ends have been updated in self.
    # iterate through the junctions.
    # check for a left overhang.
    numfuz2left = 0
    numselfleft = 0
    if not self.fuzzy_junctions[0].overlaps(fuz2.fuzzy_junctions[0],self.junction_tolerance):
      # see if we need to add sequences from fuz2
      if mode(fuz2.fuzzy_junctions[0].left.get_payload()['junc']) < mode(self.fuzzy_junctions[0].left.get_payload()['junc']):
         #print 'left over2'
         i = 0
         while not self.fuzzy_junctions[0].overlaps(fuz2.fuzzy_junctions[i],self.junction_tolerance) and i < len(fuz2.fuzzy_junctions):
           i+=1
         numfuz2left = i # number to push on from the fuz2 and increment in
         #print numfuz2left
      elif mode(fuz2.fuzzy_junctions[0].left.get_payload()['junc']) > mode(self.fuzzy_junctions[0].left.get_payload()['junc']):
         #print 'left over1'
         i = 0
         while not self.fuzzy_junctions[i].overlaps(fuz2.fuzzy_junctions[0],self.junction_tolerance) and i < len(self.fuzzy_junctions):
           i+=1
         numselfleft = i # number to increment in from self
         #print numselfleft
      else:
        sys.stderr.write("WARNING: strange case \n")
        return False
    # next we can check how long we have a run of the same 
    ind1 = numselfleft
    ind2 = numfuz2left
    overlap_size = 0
    while ind1 < len(self.fuzzy_junctions) and ind2 < len(fuz2.fuzzy_junctions) \
      and self.fuzzy_junctions[ind1].overlaps(fuz2.fuzzy_junctions[ind2],self.junction_tolerance):
      overlap_size += 1
      ind1 += 1
      ind2 += 1
    #print 'overlap size '+str(overlap_size)
    numselfright = len(self.fuzzy_junctions) - overlap_size - numselfleft
    numfuz2right = len(fuz2.fuzzy_junctions) - overlap_size - numfuz2left
    if min(numselfright,numfuz2right) != 0: 
      sys.stderr.write("WARNING: expected one of them to be zero\n")
      return False
    if min(numselfleft,numfuz2left) != 0: 
      sys.stderr.write("WARNING: expected one of them to be zero\n")
      return False
    #print numselfright
    #print numfuz2right
    #print self.fuzzy_junctions[numselfleft].overlaps(fuz2.fuzzy_junctions[numfuz2left],self.junction_tolerance)
    #print 'add'
    #Now we have what we need to go through and do some updating
    #Lets just make new fuzzy junctions
    newjuncs = []
    for i in range(0,numfuz2left):
      newjuncs.append(fuz2.fuzzy_junctions[i])
    for i in range(0,numselfleft):
      newjuncs.append(self.fuzzy_junctions[i])
    #Now we do both down the center
    range1 = range(numselfleft,overlap_size+numselfleft)
    range2 = range(numfuz2left,overlap_size+numfuz2left)    
    for i in range(0,len(range1)):
      newjuncs.append(self.fuzzy_junctions[range1[i]])
      newjuncs[-1].add_fuzzy_junction(fuz2.fuzzy_junctions[range2[i]])
      #print i
    #Make the right size
    for i in range(overlap_size+numfuz2left,overlap_size+numfuz2left+numfuz2right):
      newjuncs.append(fuz2.fuzzy_junctions[i])
    for i in range(overlap_size+numselfleft,overlap_size+numselfleft+numselfright):
      newjuncs.append(self.fuzzy_junctions[i])
    self.fuzzy_junctions = newjuncs
    #print 'adding gpd '+str(len(fuz2.gpds))+' entries'
    for g in fuz2.gpds: self.gpds.append(g)
    #print 'new entry'
    #print self.get_info_string()
    return True

  def do_add_single_exon_fuzzy_gpd(self,fuz2):
    if not self.do_add_single_exon:
      return False  # make sure we are allowed to be doing this
    #build the bounds from the average start and end
    s1 = mean(self.start.get_payload())
    e1 = mean(self.end.get_payload())
    s2 = mean(fuz2.start.get_payload())
    e2 = mean(fuz2.end.get_payload())
    l1 = e1-s1+1
    l2 = e2-s2+1
    if l1 < self.single_exon_minimum_length:
      return False
    if l2 < self.single_exon_minimum_length:
      return False
    if l1 < 1 or l2 < 1: return False #shouldn't happen
    chr1 = self.start.chr
    chr2 = self.end.chr
    if chr1 != chr2: return False #shouldn't happen
    r1 = Bed(chr1,s1-1,e1,self.dir)
    r2 = Bed(chr2,s2-1,e2,self.dir)
    over = r1.overlap_size(r2)
    if over < self.single_exon_minimum_overlap_bases:
      return False
    #print r1.get_range_string()
    #print r2.get_range_string()
    cov = min(float(over)/float(l1),float(over)/float(l2))
    if cov < self.single_exon_minimum_overlap_fraction:
      return False
    if abs(e1-e2) > self.single_exon_maximum_endpoint_distance:
      return False
    if abs(s1-s2) > self.single_exon_maximum_endpoint_distance:
      return False
    #If we're still here, we can add result
    newstart = self.start.merge(fuz2.start)
    newstart.set_payload([])
    for s in self.start.get_payload():
      newstart.get_payload().append(s)
    for s in fuz2.start.get_payload():
      newstart.get_payload().append(s)
    newend = self.end.merge(fuz2.end)
    newend.set_payload([])
    for e in self.end.get_payload():
      newend.get_payload().append(e)
    for e in fuz2.end.get_payload():
      newend.get_payload().append(e)
    self.start = newstart
    self.end = newend
    for gpd in fuz2.gpds: self.gpds.append(gpd)
    return True

  #Return true if these genepreds can be added together
  def compatible_overlap(self,fingpd):
    f1 = self
    f2 = fingpd
    #### Forget about trying zero exon cases for now
    if len(f1.fuzzy_junctions)==0 or len(f2.fuzzy_junctions)==0:
      return False
    matches = []
    for i in range(0,len(f1.fuzzy_junctions)):
      for j in range(0,len(f2.fuzzy_junctions)):
        if f1.fuzzy_junctions[i].overlaps(f2.fuzzy_junctions[j],self.junction_tolerance):
          matches.append([i,j])
    # This is our matched junctions in f1 and f2
    if len(matches)==0: return False
    # This is the number of extra exons it would take in the middle of the run (shifts)
    if len(set([x[0]-x[1] for x in matches])) != 1:  return False
    # Lets make sure all our exons are consecutive
    if len(matches) > 1:
      consec1 = list(set([matches[i+1][0]-matches[i][0] for i in range(0,len(matches)-1)]))
      consec2 = list(set([matches[i+1][1]-matches[i][1] for i in range(0,len(matches)-1)]))
      if len(consec1) != 1: return False
      if len(consec2) != 1: return False
      if consec1[0] != 1: return False
      if consec2[0] != 1: return False
    # one of them should be zero
    if not(matches[0][1] == 0 or matches[0][0] == 0):
      return False
    # and one of our last matches should be the last junction
    if not (len(f1.fuzzy_junctions)-1==matches[-1][0] or len(f2.fuzzy_junctions)-1==matches[-1][1]):
      return False
    ## because of how things are ordered our offset should always be positive
    #offset = [x[0]-x[1] for x in matches][0]
    #print 'offset: '+str(offset)
    #print matches
    #print [x[0]-x[1] for x in matches]
    #print len(set([x[0]-x[1] for x in matches]))
    #if len(matches) > 1:
    #  consec = list(set([matches[i+1][0]-matches[i][0] for i in range(0,len(matches)-1)]))
    #  print consec
    #  consec = list(set([matches[i+1][1]-matches[i][1] for i in range(0,len(matches)-1)]))
    #  print consec
    #print '---'
    return True

  def read_first(self,ingpd):
      self.gpds.append(ingpd)
      if self.use_dir: self.dir = ingpd.value('strand')
      # add fuzzy junctions
      chr = ingpd.value('chrom')
      for i in range(0,len(ingpd.value('exonStarts'))-1):
        self.fuzzy_junctions.append(FuzzyJunction(chr,ingpd.value('exonEnds')[i],ingpd.value('exonStarts')[i+1],self.dir))
      if len(ingpd.value('exonStarts')) > 1:
        self.fuzzy_junctions[0].left.get_payload()['start'] = Bed(chr,ingpd.value('txStart'),ingpd.value('txStart')+1,self.dir)
        self.fuzzy_junctions[0].left.get_payload()['start'].set_payload([])
        self.fuzzy_junctions[0].left.get_payload()['start'].get_payload().append(ingpd.value('txStart')+1)
        self.fuzzy_junctions[-1].right.get_payload()['end'] = Bed(chr,ingpd.value('txEnd')-1,ingpd.value('txEnd'),self.dir)
        self.fuzzy_junctions[-1].right.get_payload()['end'].set_payload([])
        self.fuzzy_junctions[-1].right.get_payload()['end'].get_payload().append(ingpd.value('txEnd'))
      # add fuzzy starts
      self.start = Bed(ingpd.value('chrom'),ingpd.value('txStart'),ingpd.value('txStart')+1,self.dir)
      self.start.set_payload([])
      self.start.get_payload().append(ingpd.value('txStart')+1)
      self.end = Bed(ingpd.value('chrom'),ingpd.value('txEnd')-1,ingpd.value('txEnd'),self.dir)
      self.end.set_payload([])
      self.end.get_payload().append(ingpd.value('txEnd')+1)
      # Have finished reading in the first case

class FuzzyJunction:
  # Pre: inleft is 1-indexed last exonic base on the left
  #      inright is 1-indexed first exonic base on the right
  #      direction doesn't need to be used
  def __init__(self,inchr=None,inleft=None,inright=None,indir=None):
    self.chr = inchr
    self.left = None  #range with payloads being the actual left and rights
    self.right = None
    self.dir = indir
    if inchr and inleft and inright:
      self.add_junction(inchr,inleft,inright,indir)

  # return chr, and the left and right mode as an array
  def get_mode(self):
    m1 = mode(self.left.get_payload()['junc'])
    m2 = mode(self.right.get_payload()['junc'])
    return [Bed(self.chr,m1-1,m1,self.dir),Bed(self.chr,m2-2,m2,self.dir)]
  # Find the mode of the junction and see if this overlaps
  def overlaps(self,fjun2,juntol):
    m1 = self.get_mode()
    m2 = fjun2.get_mode()
    if m1[0].direction != m2[0].direction: return False # usually they are both off
    if not m1[0].overlaps_with_padding(m2[0],juntol): return False
    if not m1[1].overlaps_with_padding(m2[1],juntol): return False
    return True

  #Right now assumes these are overlap verified prior to calling
  def add_junction(self,inchr,inleft,inright,indir=None):
    if not self.left: # this is our first one
      t1 = {}
      t1['junc'] = []
      t1['start'] = None
      self.left = Bed(inchr,inleft-1,inleft,indir)
      self.left.set_payload(t1)
      self.left.get_payload()['junc'].append(inleft)
      self.right = Bed(inchr,inright-1,inright,indir)
      t2 = {}
      t2['junc'] = []
      t2['end'] = None
      self.right.set_payload(t2)
      self.right.get_payload()['junc'].append(inright)
      return
    #Lets add this one to our current one
    if inchar != self.chr:
      sys.stderr.write("ERROR: need to be overlap checked ahead of this\n")
      sys.exit()
    if self.dir != indir:
      sys.stderr.write("ERROR: need to be overlap checked ahead of this\n")
      sys.exit()
    newfuz = FuzzyJunction(inchar,inleft,inright,indir)
    self.add_fuzzy_junction(newfuz)
  def add_fuzzy_junction(self,newfuz):
    #print 'add fuzzy'
    mergeleft = self.left.merge(newfuz.left)
    mergeleft.set_payload(self.left.get_payload())
    mergeright = self.right.merge(newfuz.right)
    mergeright.set_payload(self.right.get_payload())
    for j1 in newfuz.left.get_payload()['junc']:
      mergeleft.get_payload()['junc'].append(j1)
    for j2 in newfuz.right.get_payload()['junc']:
      mergeright.get_payload()['junc'].append(j2)
    #fix the starts
    if newfuz.left.get_payload()['start'] and not self.left.get_payload()['start']:
      mergeleft.get_payload()['start'] = newfuz.left.get_payload()['start']
    elif newfuz.left.get_payload()['start'] and self.left.get_payload()['start']:
      newrange = self.left.get_payload()['start'].merge(newfuz.left.get_payload()['start'])
      newrange.set_payload([])
      for s in self.left.get_payload()['start'].get_payload(): newrange.get_payload().append(s)
      for s in newfuz.left.get_payload()['start'].get_payload(): newrange.get_payload().append(s)
      mergeleft.get_payload()['start'] = newrange
      #print 'update left starts'
    #fix the ends
    if newfuz.right.get_payload()['end'] and not self.right.get_payload()['end']:
      mergeright.get_payload()['end'] = newfuz.right.get_payload()['end']
    elif newfuz.right.get_payload()['end'] and self.right.get_payload()['end']:
      newrange = newfuz.right.get_payload()['end'].merge(self.right.get_payload()['end'])
      newrange.set_payload([])
      for s in self.right.get_payload()['end'].get_payload(): newrange.get_payload().append(s)
      for s in newfuz.right.get_payload()['end'].get_payload(): newrange.get_payload().append(s)
      mergeright.get_payload()['end'] = newrange
      #print 'update right ends'
    # We finished the changes
    self.left = mergeleft
    self.right = mergeright
    #print 'done fuzzy'
    #print mergeleft.get_range_string()
    
class FuzzyGenePredSeparator:
  def __init__(self):
    self.junction_tolerance = 10 #bp tolerance for junction coordinate match
    self.use_dir = False # do we use direction
    self.full_length = False # if False we meld in genepreds that are subsets of others
    #self.gpds = []  #member genepreds
    return
  def set_junction_tolerance(self,juntol):
    self.junction_tolerance = juntol
  def set_use_direction(self,usedir):
    self.use_dir = usedir
  def set_full_length(self,full_length):
    self.full_length = full_length
  def get_fuzzies(self,gpds):
    outs = []
    for gpd in gpds:
      outs.append(FuzzyGenePred(gpd,use_dir=self.use_dir,full_length=self.full_length,jun_tol=self.junction_tolerance))
    return outs

def split_genepreds_by_overlap(gpds,use_dir=False):
  if len(gpds) == 1:  return [gpds]
  sets = []
  for g in gpds:  sets.append([g])
  olen = 0
  while olen != len(sets):
    olen = len(sets)
    sets = combine_down_gpd(sets,use_dir)
  #print len(sets)
  #for set1 in sets:
  #  print '  '+str(len(set1))+"\t"+set1[0].get_bed().get_range_string()
  return sets

#slow implementation
def combine_down_gpd2(sets,use_dir=False):
  for i in range(0,len(sets)):
    for j in range(i+1,len(sets)):
      if use_dir and sets[i][0].value('strand') != sets[j][0].value('strand'): continue
      match = False
      for k in sets[i]:
        for l in sets[j]:
          for m in k.range_set.get_range_list():
            for n in l.range_set.get_range_list():
              if m.overlaps(n):
                match = True
                if match: break
            if match: break
          if match: break
        if match: break
      if match: #make the change and get out
        for x in sets[j]: sets[i].append(x)
        del sets[j]
        return 1
  return 0

def combine_down_gpd(sets,use_dir=False):
  osets = []
  while len(sets) > 0:
    curr = sets.pop(0)
    # see if curr over laps sets
    match = False
    for i in range(0,len(osets)):
      if use_dir and curr[0].value('strand') != osets[i][0].value('strand'): continue
      for m in curr:
        for n in osets[i]:
          if m.overlaps(n): match = True
          break
        if match: break
      if match:
        for m in curr: osets[i].append(m)
        break
    if not match:
      osets.append(curr)
  return osets
      
def mode(list):
  return max(set(list),key=list.count)

# Pre: a list of genepreds and a junction tolerance
# Post:  a list fuzzy genepreds where they have been combined where appropriate
def gpd_list_to_combined_fuzzy_list(gpds,juntol,full_length=False,use_dir=False):
  fgs = FuzzyGenePredSeparator()
  fgs.set_junction_tolerance(juntol)
  fgs.set_use_direction(use_dir)
  fgs.set_full_length(full_length)
  splitgpds = split_genepreds_by_overlap(gpds,use_dir=False)
  results = []
  cnt = 0
  for gset in splitgpds:
    cnt += len(gset)
    fzs = fgs.get_fuzzies(gset)
    #see if we can add these
    prevlen = 0
    while len(fzs) != prevlen:
      prevlen = len(fzs)
      fzs = combine_down_fuzzies(fzs)
    for o in fzs:  results.append(o)
  #print 'worked on '+str(cnt)+' genepreds'
  return results

def combine_down_fuzzies(fzs):
  if len(fzs) == 1:  return fzs
  outs = []
  while len(fzs) > 0:
    curr = fzs.pop(0)
    combined = False
    for i in range(0,len(outs)):
      combined = outs[i].add_fuzzy_gpd(curr)
      if combined: break
    if not combined:
      outs.append(curr)
  #for fz in outs:
  #  print fz.get_info_string()
  return outs

def mean(list):
  if len(list) == 1: return list[0]
  return int(float(sum(list))/float(len(list)))
