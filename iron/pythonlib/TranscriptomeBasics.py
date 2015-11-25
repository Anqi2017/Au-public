import sys, re, hashlib, json, random
import GenePredBasics, SequenceBasics

# Transcriptome is a set of genepred entries
# with the corresponding fasta file.
# alternatively, you can read in a serialized transcriptome.
#
# You can further define a transcriptome file with an expression file
# This file can be of the form of a TSV
#
class Transcriptome:
  def __init__(self):
    #self.transcript_names = {}
    self.transcripts = {}
    #self.gpds = {}
    self.expression = None
    self.ref_hash = None
  def set_reference_genome_dictionary(self,indict):
    self.ref_hash = indict
    return
  def add_expression(self,inname,exp):
    if not self.expression:
      self.expression = IsoformExpression()
      for name in self.transcripts: self.expression.add_expression(name,0)
    self.expression.add_expression(inname,exp)

  def get_serialized(self):
    jexpress = None
    if self.expression is not None:
      jexpress = self.expression.get_serialized()
    #return json.loads(
    todump = [self.transcript_names,self.transcripts,self.gpds,jexpress]
    return json.dumps(todump)
  def read_serialized(self,input):
    [self.transcript_names,self.transcripts,self.gpds,jexpress] = json.loads(input)
    if jexpress:
      self.expression = IsoformExpression()
      self.expression.read_serialized(jexpress)
    else:
      self.expression = None
  def add_genepred_line(self,inline):
    if not self.ref_hash:  
      sys.stderr.write("ERROR: Must assign a reference genome dictionary first\n")
      sys.exit()
    gpd = GenePredBasics.GenePredEntry(inline)
    if gpd.value('name') in self.transcripts:
      sys.stderr.write("WARNING: "+inline+" transcript was already set\n")
    seq = ''
    for i in range(0,gpd.value('exonCount')):
      seq += self.ref_hash[gpd.value('chrom')][gpd.value('exonStarts')[i]:gpd.value('exonEnds')[i]].upper()
    if gpd.value('strand') == '-': seq = SequenceBasics.rc(seq)
    self.transcripts[gpd.value('name')] = seq
    return    

  # This is depreciated
  #def read_from_fasta_and_genepred(self,genomefastafile,genepredfile):
  #  # read in our genome
  #  seen_names = {}
  #  seen_coords = {}
  #  genepred = {}
  #  with open(genepredfile) as inf:
  #    for line in inf:
  #      if re.match('^#',line): continue
  #      e = GenePredBasics.line_to_entry(line)
  #      hexcoord = hashlib.sha1(e['chrom']+"\t"+e['strand'] + "\t" + str(e['exonStarts'])+"\t" + str(e['exonEnds'])).hexdigest()
  #      dupname = 0
  #      dupcoord = 0
  #      if hexcoord in seen_coords:
  #        sys.stderr.write("Warning "+ e['name'] + " " + e['gene_name'] + " exists at identical coordinates as another entry\n")
  #        dupcoord = 1
  #      seen_coords[hexcoord] = 1
  #      currname = e['name']
  #      if e['name'] in seen_names:
  #        if dupcoord == 1:
  #          sys.stderr.write("skipping perfect duplicate of "+e['name']+"\n")
  #          continue
  #        newname = e['name'] + "."+str(len(seen_names[e['name']])+1)
  #        currname = newname
  #        seen_names[e['name']].append(newname)
  #        sys.stderr.write("Warning "+ e['name'] + " " + e['gene_name'] + " is a duplicate name.. renaming to "+newname+ "\n")
  #        dupname = 1
  #      else:
  #        seen_names[e['name']] = []
  #        seen_names[e['name']].append(e['name'])
  #      genepred[currname] = e
  #
  #  #print "reading names and locs"             
  #  ref = SequenceBasics.read_fasta_into_hash(genomefastafile)
  #  #print "converting sequences"
  #  for transcript in genepred:
  #    e = genepred[transcript]
  #    if e['chrom'] in ref:
  #      seq = ''
  #      self.transcript_names[transcript] = genepred[transcript]['name']
  #      for i in range(0,e['exonCount']):
  #        seq += ref[e['chrom']][e['exonStarts'][i]:e['exonEnds'][i]]
  #      if e['strand'] == '-': seq = SequenceBasics.rc(seq)
  #      self.transcripts[transcript] = seq.upper()
  #      self.gpds[transcript] = e

  # Pre: Expression must have been set
  # Post: Returns a random transcript name
  def get_random_by_expression(self):
    return self.expression.get_random_by_expression()    
  def get_uniform_random(self):
    tnames = self.transcripts.keys()
    tnum = len(tnames)
    rnum = random.randint(0,tnum-1)
    return tnames[rnum]
  # Default to random by expression if its set
  def get_random(self):
    if self.expression: return self.get_random_by_expression()
    return self.get_uniform_random()
  def get_sequence(self,name):
    if name not in self.transcripts:
      sys.stderr.write("ERROR: "+name+" not in transcripts\n")
      sys.exit()
    return self.transcripts('name')


# Class holds the isoform names and expression values
# And also has functions for randomly getting an isoform name
# either by uniform distribution or 
class IsoformExpression:
  def __init__(self):
    self.expression = {}
    self.total_expression = None
    self.names = None
    return
  # Pre:  TSV with <transcript name> <expression level>
  def read_tsv(self,filename):
    with open(filename) as inf:
      for line in inf:
        f = line.rstrip().split("\t")
        self.expression[f[0]]=float(f[1])
    self.update_expression()
  def get_expression(self,transcript_name):
    if transcript_name not in self.expression: 
      sys.stderr.write("ERROR: "+transcript_name+" not in expression")
      sys.exit()
    return self.expression[transcript_name]
  # Add a single expression value
  def add_expression(self,transcript_name,expression):
    self.expression[transcript_name] = expression
    self.update_expression()
  def read_serialized(self,instring):
    [self.expression] = json.loads(instring)
    self.get_total_expression()
  def get_serialized(self):
    return json.dumps([self.expression])
  def get_random_by_expression(self):
    rnum = random.random()
    total = 0
    for name in self.names:
      total += self.expression[name]/self.total_expression
      if rnum < total:
        return name
    return name
  def get_uniform_random(self):
    rnum = random.randint(0,len(self.names)-1)
    return self.names[rnum]

  def update_expression(self):
    self.names = sorted(self.expression.keys())
    self.total_expression = sum([self.expression[x] for x in self.expression])
    
