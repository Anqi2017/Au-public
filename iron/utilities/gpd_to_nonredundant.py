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
  downsampcount = 0
  if args.threads > 1:
    p = Pool(processes=args.threads)
  while True:
    lcount += 1
    locus = gpdls.read_locus()
    if not locus: break
    garuntee = 500
    if len(locus) > args.downsample+garuntee:
      downsampcount += 1
      sys.stderr.write("downsampling "+str(downsampcount)+"\n")
      nlocus = sorted(locus,key=lambda x: x.length(), reverse=True)
      locus = []
      for x in nlocus[0:garuntee]: locus.append(x)
      mlocus = nlocus[garuntee:]
      random.shuffle(mlocus)
      for x in mlocus[0:args.downsample]: locus.append(x)
    if args.threads > 1:
      p.apply_async(process_locus,args=(lcount,locus,args),callback=do_results)
    else:
      r = process_locus(lcount,locus,args)
      do_results(r)
  if args.threads > 1:
    p.close()
    p.join()
  if args.output_original_table:
    args.output_original_table.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_results(outs):
  gpdlines,tablelines,location = outs
  if gpdlines == '': return
  global glock
  global warning_count
  global locus_count
  locus_count += 1
  glock.acquire()
  sys.stderr.write(str(locus_count)+" loci  "+location+"            \r")
  of_main.write(gpdlines)
  if of_table:
    of_table.write(tablelines)
  glock.release()

def process_locus(lcount,locus,args):
    # replace locus with a nonredundant locus set
    global warning_count
    global glock
    [nrlocus,nrlocuskey] = get_nr_locus(locus)
    nrfuzzykey = {}
    location = None
    for num in nrlocuskey:
      v = greedy_combine_down_fuzzies([FuzzyGenePred(x,juntol=args.junction_tolerance*2) for x in nrlocuskey[num]])
      if len(v) > 1:
        if args.verbose:  sys.stderr.write("WARNING expected only 1 fuzzy genepred\n")
        glock.acquire()
        warning_count += 1
        glock.release()
      nrfuzzykey[num] = v[0]
    [subset,compatible] = do_locus(lcount,nrlocus,args)
    if args.predict:
      do_prediction(compatible,args,nrfuzzykey,location)
    else:
      return do_reduction(subset,args,nrfuzzykey,location)

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
  
def do_locus(lcount,locus,args):
  fname = args.tempdir+'/'+str(lcount)+'.bed'
  cmd1 = "bedtools sort -i -"
  of = open(fname,'w')
  p1 = Popen(cmd1.split(),stdin=PIPE,stdout=of)
  for rnum in range(0,len(locus)):
    e = locus[rnum]
    #convert each genpred to junction information
    if e.entry['exonCount'] < 2: # single exon
      continue
    # have multiple exons
    #print e.entry['gene_name']
    for i in range(0,e.entry['exonCount']-1):
      jnum = i
      # now get junctions
      bstr1 = ''
      bstr1 += e.entry['chrom']+"\t"
      leftleft = e.entry['exonEnds'][i]-args.junction_tolerance-1
      if leftleft < 0: leftleft = 0
      bstr1 += str(leftleft)+"\t"
      bstr1 += str(e.entry['exonEnds'][i]+args.junction_tolerance)+"\t"
      bstr1 += str(rnum)+"\t"
      bstr1 += str(jnum)+"\t"
      bstr1 += "l"+"\t"
      bstr1 += str(e.entry['exonEnds'][i])+"\t"
      bstr1 += str(e.entry['exonCount']-1)
      bstr2 = ''
      bstr2 += e.entry['chrom']+"\t"
      rightright = e.entry['exonStarts'][i+1]-args.junction_tolerance
      if rightright < 0: rightright = 0
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
    #if intcount%10000 == 0: sys.stderr.write(str(intcount)+"        \r")
    f = line.rstrip().split("\t")
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
    if r1 == r2: continue #skip match == 
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
        if len(jres[r1][j1][r2]) != 2:
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
        iscomp = is_compatible(reads[r1][r2],jtotals[r1],jtotals[r2])
        if iscomp:
          if r1 not in cout: cout[r1] = set()
          cout[r1].add(r2) 
  return [sout,cout]

#check if r2 could be added to r1
def is_compatible(ematch,c1,c2):
  if c1 == 1 or c2 == 1: return True
  if ematch[0][0] !=0 or ematch[0][1] != 0: return False 
  if ematch[-1][0] != c1-1 or ematch[-1][1] != c2-1: return False
  maxgap1 = max([ematch[i+1][0]-ematch[i][0] for i in range(0,len(ematch)-1)])
  maxgap2 = max([ematch[i+1][1]-ematch[i][1] for i in range(0,len(ematch)-1)])
  if maxgap1 == 1 and maxgap2 == 1:  return True
  return False

#check if r2 is a subset of r1
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
  parser=argparse.ArgumentParser(description="")
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--output_original_table',help="A table to translate new names to original evidence\n")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  parser.add_argument('-j','--junction_tolerance',default=0,type=int)
  parser.add_argument('-v','--verbose',action='store_true')
  parser.add_argument('--downsample',type=int,default=2000,help="Maximum read depth at locus. sample down to random subset this size")
  parser.add_argument('--predict',action='store_true',help="build out longer reads based on the inputs")
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
