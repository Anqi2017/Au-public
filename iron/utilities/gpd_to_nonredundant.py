#!/usr/bin/python
import argparse, sys, os, random
from shutil import rmtree
from multiprocessing import cpu_count, Pool, Lock
from tempfile import mkdtemp, gettempdir
from GenePredBasics import GenePredLocusStream, GenePredEntry
from subprocess import Popen, PIPE
import GraphBasics
from GenePredFuzzyBasics import FuzzyGenePred, greedy_combine_down_fuzzies

glock = Lock()
of_main = None
of_table = None
warning_count = 0
locus_count = 0
downsampcount = 0

def main():
  #do our inputs
  args = do_inputs()
  if args.output_original_table: args.output_original_table = open(args.output_original_table,'w')
  global of_table 
  of_table = args.output_original_table
  global of_main
  of_main = args.output
  gpdls = GenePredLocusStream(args.input)
  lcount = 0
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for locus in gpdls:
    lcount += 1
    # When a locus too many reads use random downsampling
    # to reduce to a more manageable number of reads
    if len(locus) > args.downsample:
      locus = downsample(locus,args.downsample,500)
    #Execute Per Locus analysis with results passed to do_results() callback
    if args.threads > 1:
      p.apply_async(process_locus,args=(lcount,locus,args),callback=do_results)
    else:
      r = process_locus(lcount,locus,args)
      do_results(r)
  if args.threads > 1:
    p.close()
    p.join()
  #Close the table output if we are making it
  if args.output_original_table:
    args.output_original_table.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

# Pre: locus - List of gpds in same locus
#      downsize - new size of output
#      guarantee - keep the X longest alignments then also
#                  the rest will be randomly sampled down
# Post: locus array with 'guarantee' longest alignments and
#       'downsize' - 'guarantee' randomly sampled alignments 
# Modifies: Warns STDERR that downsampling is occuring
#           downsampcount global is updated
def downsample(locus,downsize,guarantee):
  global downsampcount
  global glock
  glock.acquire()
  downsampcount += 1
  glock.release()
  if downsampcount < 100:
    sys.stderr.write("\ndownsampling "+str(downsampcount)+" " +str(len(locus))+"\n")
  nlocus = sorted(locus,key=lambda x: x.length(), reverse=True)
  locus = []
  for x in nlocus[0:guarantee]: locus.append(x)
  mlocus = nlocus[guarantee:]
  random.shuffle(mlocus)
  for x in mlocus[0:downsize]: locus.append(x)
  return locus

# Callback to output results
# Pre: array contains gpdlines, tablelines, location
#      gpdlines is the gpd output
#      tablelines is table otuput
#      location is just to help keep track of progress
# Post: Writes to global outputs of_main and of_table
# Modifies: Writes to outputs defined in globals of_main and of_table
#           Iterates global locus_count
#           Uses global glock Lock()
#           Status update to STDERR
def do_results(outs):
  if not outs: return
  gpdlines,tablelines,location = outs
  if gpdlines == '': return
  global glock
  global warning_count
  global locus_count
  global of_main
  global of_table
  locus_count += 1
  glock.acquire()
  sys.stderr.write(str(locus_count)+" loci  "+location+"            \r")
  of_main.write(gpdlines)
  if of_table:
    of_table.write(tablelines)
  glock.release()

# Run the processing of the locus
#   Executes either a 'do_reduction' or 'do_prediction'
#   processing of data based on args
# Pre: lcount - the count for this locus
#      locus - array of GPDs
#      args - argparse args
# Post: Returns a list for use in the callback
#       0 - locus range
#       1 - combined gpd lines
#       2 - table of new gpd names to combined gpd line names
def process_locus(lcount,locus,args):
  # replace locus with a nonredundant locus set
  #Simplify data by reducing exact duplicates
  [nrlocus,nrlocuskey] = get_nr_locus(locus)
  #Further simplify data by putting like gpds into fuzzy gpds
  nrfuzzykey = get_nr_fuzzy(nrlocuskey,args)
  location = None
  [subset,compatible] = do_locus(lcount,nrlocus,args)
  if args.predict:
    return do_prediction(compatible,args,nrfuzzykey,location)
  else:
    return do_reduction(subset,args,nrfuzzykey,location)

# Simplify same GPDS by making them fuzzy GPDs
# Pre: nrlocuskey has a dict and the args
#      The dict contains sets of like gpds
# Post: after having combined down in fuzzy gpds return them in a dict
def get_nr_fuzzy(nrlocuskey,args):
  global warning_count
  global glock
  nrfuzzykey = {}
  for num in nrlocuskey:
    # Create FuzzyGenePreds out of all of the GPDs
    # And reduce the sets of many gpds down to just a single fuzzy gpd representing each specific gpd.
    v = greedy_combine_down_fuzzies([FuzzyGenePred(x,juntol=args.junction_tolerance*2) for x in nrlocuskey[num]])
    if len(v) > 1:
      if args.verbose:  sys.stderr.write("WARNING expected only 1 fuzzy genepred\n")
      glock.acquire()
      warning_count += 1
      glock.release()
    nrfuzzykey[num] = v[0]
  return nrfuzzykey

# Build up a set of predicted gpds
# Pre:
#   compatible - hash to say if one read index is compatible with another
#              compatible[r1][r2] r1 is compatible with r2
#   args - argparse inut
#   nrfuzzykey - by index each genepred stored in fuzzy format
#   location - Not sure if this is necessary
# Post:
def do_prediction(compatible,args,nrfuzzykey,location):
    #if len(compatible.keys()) == 0: return None
    #all reads could be standing alone version
    families = []
    for num in nrfuzzykey:
      families.append(nrfuzzykey[num])
      nrfuzzykey[num].params['proper_set'] = False #partial overlap is enough
    #get_compatible_evidence(compatible,nrfuzzykey,args)
    for i in compatible:
      for j in compatible[i]:
        #see if its already in there
        g1lines = set()
        for g1 in nrfuzzykey[i].gpds: g1lines.add(g1.get_line())
        repeat = False
        for g2 in nrfuzzykey[j].gpds:
          if g2.get_line() in g1lines:
            repeat = True
            break
        if not repeat: continue
        together = nrfuzzykey[i].concat_fuzzy_gpd(nrfuzzykey[j])
        if together:
          families.append(together)
    # now we need to find any duplicate entries and combine them
    newfam = []
    beforefam = len(families)
    while len(families) > 0:
      fam = families.pop(0)
      remaining = []
      for i in range(0,len(families)):
        if fam.is_equal_fuzzy(families[i]):
          added = fam.add_fuzzy_gpd(families[i])
          if not added:
            sys.stderr.write("WARNING NOT SURE WHY NOT ADDED EQUAL\n")
          fam = added
        else: remaining.append(families[i])
      families = remaining
      newfam.append(fam)
    families = newfam
    afterfam = len(families)

    # Replace the family with a set where we haven't used the same gpd line twice
    # This may damage the fuzzy object
    for i in range(0,len(families)):
      gset = set()
      for g in families[i].gpds:  
        gset.add(g.get_line())
      families[i].gpds  = [GenePredEntry(x) for x in gset]
    #  sys.stderr.write("\n\ncahnged from "+str(beforefam)+"\t"+str(afterfam)+"\n\n")
    gpdlines = ""
    tablelines = ""
    # find gpds not in the graph... 
    for fz in families:
      info = fz.get_info_string()
      gpdline = fz.get_genepred_line()
      #print '&&&&&&&&&&&&&&&&'
      #print gpdline
      #print fz.get_info_string()
      #print '&&&&&&&&&&&&&&&&'
      gpd = GenePredEntry(gpdline)
      if not gpd.is_valid(): 
        sys.stderr.write("WARNING: invalid genepred entry generated\n"+gpdline+"\n"+fz.get_info_string()+"\n")
        gpd = sorted(fz.gpds, key=lambda x: x.get_exon_count(), reverse=True)[0] #just grab one that has all the exons
        fz = FuzzyGenePred(gpd,juntol=args.junction_tolerance*2)
        gpdline = fz.get_genepred_line()
        if not gpd.is_valid():
          sys.stderr.write("WARNING: still problem skilling\n")
          continue
      gpdlines += gpdline+"\n"
      if args.output_original_table:
        name = gpd.entry['name']
        for g in fz.gpds:
          tablelines+=name+"\t"+g.entry['name']+"\n"
      grng = gpd.get_bed()
      grng.direction = None
      if not location: 
        location = grng
      location = location.merge(grng)
    locstring = ''
    if location:  locstring = location.get_range_string()
    return [gpdlines, tablelines, locstring]

def do_reduction(subset,args,nrfuzzykey,location):
    seen = set()
    for i in subset:
      seen.add(i)
      for j in subset[i]:  seen.add(j)
    singles = []
    for num in nrfuzzykey:
      if num not in seen:
        singles.append(num)
    #if len(subset.keys()) == 0 and len(compatible.keys()) == 0: return
    families = get_subset_evidence(subset,nrfuzzykey,args)
    gpdlines = ""
    tablelines = ""
    for num in singles:
      families.append(nrfuzzykey[num])
    # find gpds not in the graph... 
    for fz in families:
      info = fz.get_info_string()
      gpdline = fz.get_genepred_line()
      #print '&&&&&&&&&&&&&&&&'
      #print gpdline
      #print fz.get_info_string()
      #print '&&&&&&&&&&&&&&&&'
      gpd = GenePredEntry(gpdline)
      if not gpd.is_valid(): 
        sys.stderr.write("WARNING: invalid genepred entry generated\n"+gpdline+"\n"+fz.get_info_string()+"\n")
        gpd = sorted(fz.gpds, key=lambda x: x.get_exon_count(), reverse=True)[0] #just grab one that has all the exons
        fz = FuzzyGenePred(gpd,juntol=args.junction_tolerance*2)
        gpdline = fz.get_genepred_line()
        if not gpd.is_valid():
          sys.stderr.write("WARNING: still problem skilling\n")
          continue
      gpdlines += gpdline+"\n"
      if args.output_original_table:
        name = gpd.entry['name']
        for g in fz.gpds:
          tablelines+=name+"\t"+g.entry['name']+"\n"
      grng = gpd.get_bed()
      grng.direction = None
      if not location: 
        location = grng
      location = location.merge(grng)
    locstring = ''
    if location:  locstring = location.get_range_string()
    return [gpdlines, tablelines, locstring]

# Takes a directional graph subset
# Takes a nrfuzzykey that is the fuzzygenepred of everything keyed in the graph
# Makes an actual graph and traverses it in order to combine values
def get_subset_evidence(subset,nrfuzzykey,args):
    #if len(subset.keys()) == 0: continue
    # get the read numbers
    global warning_count
    global glock
    rnames = set()
    for i in subset:
      rnames.add(i)
      for j in subset[i]: rnames.add(j)
    g = GraphBasics.Graph()
    # create nodes
    rnodes = {}
    for i in rnames:  
      rnodes[i]=GraphBasics.Node([i])
      g.add_node(rnodes[i])
    # add edges
    for i in subset:
      for j in subset[i]:
        g.add_edge(GraphBasics.Edge(rnodes[i],rnodes[j]))
    # done with graph!
    g.merge_cycles()
    #print "============"
    #print g.get_status_string().rstrip()
    roots = g.get_roots()
    families = []
    for i in range(0,len(roots)):  
      kids = g.get_children(roots[i])
      family = []
      for p in roots[i].get_payload(): family.append(p)
      for kid in kids:
        for p in kid.get_payload(): family.append(p)
      base = nrfuzzykey[family[0]]
      if len(family) > 0:
        for id in family[1:]:
          fz = nrfuzzykey[id]
          updated = base.add_fuzzy_gpd(fz)
          if not updated:
            if args.verbose:  sys.stderr.write("WARNING: not able to add child as expected\n")
            glock.acquire()
            warning_count+=1
            glock.release()
          else:
            base = updated # make a new base to add children to
      families.append(base)
    return families


def make_sets_from_subsets(subset):
  #start with mutual subsets (or equals)
  mutual = {}
  rnames  = set()
  for r1 in subset:
    rnames.add(r1)
    for r2 in subset[r1]:
      rnames.add(r2)
      if r2 in subset:
        if r1 in subset[r2]:
          #its mutual subset
          if r1 not in mutual: mutual[r1] = set()
          mutual[r1].add(r2)
          if r2 not in mutual: mutual[r2] = set()
          mutual[r2].add(r1)
  #found mutuals
  rsets = {}
  for rname in rsets:
    rsets[rname] = set()
  
  #see if theres anything without a subset
  #nomutual = set()
  for rname in rnames:
    if rname not in subset:
      print 'no mutual for '+str(rname)
  print 'found mutuals'

# returns a representative gpd and all gpds with those junctions keyed by the representative index
def get_nr_locus(locus):
  nrlocus = []
  nrlocuskey = {}
  seen = {}
  for e in locus:
    if e.entry['exonCount'] == 1: continue
    js = e.entry['chrom'] +":"+ \
         str( [e.entry['exonEnds'][i]+e.entry['exonStarts'][i+1]+1 \
           for i in range(0,len(e.entry['exonStarts'])-1)])
    if js not in seen: seen[js] = []
    seen[js].append(e)
  ind  = 0
  for js in seen:
    nrlocus.append(seen[js][0])
    nrlocuskey[ind] = seen[js]
    ind += 1
  return [nrlocus, nrlocuskey]

# This is where we make a list of compatible gpds
# or gpds that are a proper subsets
def do_locus(lcount,locus,args):
  fname = args.tempdir+'/'+str(lcount)+'.bed'
  cmd1 = "bedtools sort -i -"
  of = open(fname,'w')
  p1 = Popen(cmd1.split(),stdin=PIPE,stdout=of,bufsize=1)
  for rnum in range(0,len(locus)):
    e = locus[rnum]
    #print e
    #convert each genpred to junction information
    if e.entry['exonCount'] < 2: # single exon
      continue
    for i in range(0,e.entry['exonCount']-1):
      jnum = i
      # now get junctions
      # we make two bed entries to be sorted that have
      # structure
      # <chr> <jtolstart> <jtolend> <read-num> <junc-num> <l> <jexact-l> <total exons>
      # <chr> <jtolstart> <jtolend> <read-num> <junc-num> <r> <jexact-r> <total exons>
      bstr1 = ''
      bstr1 += e.entry['chrom']+"\t"
      leftleft = max(0,e.entry['exonEnds'][i]-args.junction_tolerance-1)
      bstr1 += str(leftleft)+"\t"
      bstr1 += str(e.entry['exonEnds'][i]+args.junction_tolerance)+"\t"
      bstr1 += str(rnum)+"\t"
      bstr1 += str(jnum)+"\t"
      bstr1 += "l"+"\t"
      bstr1 += str(e.entry['exonEnds'][i])+"\t"
      bstr1 += str(e.entry['exonCount']-1)
      bstr2 = ''
      bstr2 += e.entry['chrom']+"\t"
      rightright = max(0,e.entry['exonStarts'][i+1]-args.junction_tolerance)
      bstr2 += str(rightright)+"\t"
      bstr2 += str(e.entry['exonStarts'][i+1]+args.junction_tolerance+1)+"\t"
      bstr2 += str(rnum)+"\t"
      bstr2 += str(jnum)+"\t"
      bstr2 += "r"+"\t"
      bstr2 += str(e.entry['exonStarts'][i+1]+1)+"\t"
      bstr2 += str(e.entry['exonCount']-1)
      p1.stdin.write(bstr1+"\n")
      p1.stdin.write(bstr2+"\n")
  p1.communicate()
  of.close()
  jres = {}
  cmd2 = "bedtools intersect -a "+fname+" -b "+fname+" -wo -sorted"
  p2 = Popen(cmd2.split(),stdout=PIPE)
  jtotals = {}
  intcount = 0
  for line in p2.stdout:
    intcount += 1
    f = line.rstrip().split("\t")
    # The position doesn't matter much here only that a match occurred.
    r1 = int(f[3])
    r2 = int(f[11])
    j1 = int(f[4])
    j2 = int(f[12])
    s1 = f[5]  #side l or r
    s2 = f[13]
    t1 = int(f[7]) # total exons
    t2 = int(f[15])
    if r1 not in jtotals: jtotals[r1] = t1
    if r2 not in jtotals: jtotals[r2] = t2
    #if r1 >= r2: continue #skip match == or redundant 
    if r1 == r2: continue #skip match == of read self alignment
    if s1 != s2: continue #skip if not correct left or right 
    if r1 not in jres: 
      jres[r1] = []
      for i in range(0,t1):
        junc = {}
        jres[r1].append(junc)
    if r2 not in jres[r1][j1]:
      jres[r1][j1][r2] = {}
    jres[r1][j1][r2][s1] = j2 
  p2.communicate()
  os.remove(fname)
  # now lets examine the overlaps
  reads = {}
  for r1 in jres:
    for j1 in range(0,len(jres[r1])):
      for r2 in jres[r1][j1]:
        if len(jres[r1][j1][r2]) != 2: # must have both 'l' and 'r'
          continue
        elif jres[r1][j1][r2]['l'] != jres[r1][j1][r2]['r']:
          continue
        # if we are still in here we have a junction match
        j2 = jres[r1][j1][r2]['r']
        match = [j1,j2]
        if r1 not in reads:
          reads[r1] = {}
        if r2 not in reads[r1]:
          reads[r1][r2] = []
        reads[r1][r2].append(match)
  sout = {}
  cout = {}
  for r1 in reads:
    for r2 in reads[r1]:
        #print str(r1) + "\t" + str(r2) + "\t" + str(reads[r1][r2])
        issub = is_subset(reads[r1][r2],jtotals[r1],jtotals[r2])
        if issub: 
          #print str(r1) + "\t" + str(r2)
          if r1 not in sout:  sout[r1] = set()
          sout[r1].add(r2)
        iscomp = is_compatible(reads[r1][r2],jtotals[r1],jtotals[r2],args)
        if iscomp:
          if r1 not in cout: cout[r1] = set()
          cout[r1].add(r2) 
  return [sout,cout]

#check if r2 could be added to r1
def is_compatible(ematch,c1,c2,args):
  #args.minimum_compatible_junctions
  #print ematch
  #if c1 == 1 or c2 == 1: return True
  if len(ematch) < args.minimum_compatible_junctions: return False #must have at least the number of junctions we will call compatible
  if len(ematch) == 0: return False
  if min(ematch[0][0],ematch[0][1]) != 0: return False #one of them must have a start
  if ematch[-1][0] != c1-1 and ematch[-1][1] != c2-1: return False #one of them must have an end
  #print 'stillin'
  #print ematch
  #print c1
  #print c2
  maxgap1 = 1
  maxgap2 = 1
  if len(ematch) > 1:
    maxgap1 = max([ematch[i+1][0]-ematch[i][0] for i in range(0,len(ematch)-1)])
    maxgap2 = max([ematch[i+1][1]-ematch[i][1] for i in range(0,len(ematch)-1)])
  if maxgap1 == 1 and maxgap2 == 1:  return True
  return False

#check if r2 is a subset of r1
# ematch - match array
# c1 - Exon count for r1
# c2 - Exon count for r2
def is_subset(ematch,c1,c2):
  if c2 == 1:  # if we only have one junction in r2, its a subset of r1
    #print 'proper subset'
    return True
  if ematch[0][1] != 0:
    # r2 has more exons laying to the left of the matches
    #print 'no proper subset'
    return False
  #compare number of exons remaining
  if ematch[-1][1] != c2-1:
    return False
  #print str(ematch) + "\t" + str(c1) + "\t" + str(c2)
  maxgap1 = max([ematch[i+1][0]-ematch[i][0] for i in range(0,len(ematch)-1)])
  maxgap2 = max([ematch[i+1][1]-ematch[i][1] for i in range(0,len(ematch)-1)])
  if maxgap1 != 1 or maxgap2 != 1:
    #print 'gap troubles'
    return False
  return True

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Position sorted gpd",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--output_original_table',help="A table to translate new names to original evidence\n")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  parser.add_argument('-j','--junction_tolerance',default=0,type=int,help="how many bases to search and combine junctions")
  parser.add_argument('-v','--verbose',action='store_true')
  parser.add_argument('--downsample',type=int,default=2000,help="Maximum read depth at locus. sample down to random subset this size")
  parser.add_argument('--predict',action='store_true',help="build out longer reads based on the inputs")
  parser.add_argument('--minimum_compatible_junctions',default=2,type=int,help="require at least this many junctions overlapped to consider two gpds compatible")
  args = parser.parse_args()
  # Setup inputs 
  if args.input == '-':
    args.input = sys.stdin
  else:
    args.input = open(args.input)
  # Setup outputs
  if args.output:
    args.output = open(args.output,'w')
  else:
    args.output = sys.stdout
  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

if __name__=="__main__":
  main()
