import re, sys, copy, subprocess, hashlib
import SequenceBasics, RangeBasics
from FileBasics import GenericFileReader

# new version of genepred_basics
# This is a class for working with genepred
# it can be seralized just by converting it back to a gpd line
class GenePredEntry:
  def __init__(self,inline=None):
    self.entry = None
    self.range_set = None
    self.locus_range = None
    self.junctions = None
    if inline:
      self.line_to_entry(inline)
  def get_bed(self):
    return RangeBasics.Bed(self.entry['chrom'],self.entry['txStart'],self.entry['txEnd'],self.entry['strand'])
  def get_smoothed(self,num):
    n = GenePredEntry(entry_to_line(smooth_gaps(self.entry,num)))
    return n
  def get_exon_count(self):
    return len(self.entry['exonStarts'])
  def get_line(self):
    return entry_to_line(self.entry)
  def length(self):
    length = 0
    for i in range(0,len(self.entry['exonStarts'])):
      length += self.entry['exonEnds'][i] - self.entry['exonStarts'][i]
    return length
  def line_to_entry(self,line):
    self.entry = line_to_entry(line)
    self.calculate_range_set()
    ecount = len(self.entry['exonStarts'])
    self.locus_range = RangeBasics.GenomicRange(self.entry['chrom'],\
                                                self.entry['exonStarts'][0]+1,\
                                                self.entry['exonEnds'][ecount-1])
    self.calculate_junctions()

  # make a range dictionary for each of these
  def calculate_range_set(self,use_dir=False):
    grd = RangeBasics.GenomicRangeDictionary()
    e = self.entry
    dir = None
    if use_dir:  dir = e['strand']
    for j in range(0,len(e['exonStarts'])):
      gr = RangeBasics.Bed(e['chrom'],e['exonStarts'][j],e['exonEnds'][j],dir)
      gr.set_payload(j)
      grd.add(gr) #put the exon number as the payload of the range dictionary because why not be able to keep them in order if we ever want to.
    self.range_set = grd
  def get_first_exon_genomic_range(self):
    for m in self.range_set.members:
      if self.entry['strand'] == '+' and m.get_payload() == 0:
        return m
      elif self.entry['strand'] == '-' and m.get_payload() == len(self.entry['exonStarts'])-1:
        return m
    sys.stderr.write("problem finding a start\n")
    return None
  def get_last_exon_genomic_range(self):
    for m in self.range_set.members:
      if self.entry['strand'] == '-' and m.get_payload() == 0:
        return m
      elif self.entry['strand'] == '+' and m.get_payload() == len(self.entry['exonStarts'])-1:
        return m
    sys.stderr.write("problem finding a start\n")
    return None
  def calculate_junctions(self):
    alljun = []
    for i in range(0,len(self.entry['exonStarts'])-1):
      jun=str(self.entry['chrom'])+':'+str(self.entry['exonEnds'][i])+','+str(self.entry['chrom'])+':'+str(self.entry['exonStarts'][i+1]+1)
      alljun.append(jun)
    self.junctions = alljun
  def value(self,vname):
      return self.entry[vname]
  #pre: another genepred
  #post: true if overlaps false if not
  def overlaps(self,gpd2,use_dir=False):
    if not self.range_set: self.calculate_range_set()
    if not gpd2.range_set: gpd2.calculate_range_set()
    #see if they overlap
    if use_dir and self.value('strand') != gpd2.value('strand'): return False
    if self.value('chrom') != gpd2.value('chrom'): return False
    for m in self.range_set.get_range_list():
      for n in gpd2.range_set.get_range_list():
        if m.overlaps(n): return True
    return False
# Compare two genepred enetry classes
# Requires a 1 to 1 mapping of exons, so if one exon overlaps two of another
# it will be a mismatch.
#   [first exon required overlap, internal exon required overlap, last exon required overlap]
#
class GenePredComparison:
  def __init__(self):
    self.overlap_requirement = [1,1,1] # require perfect first, internal, and last exon overlap by default
    self.output = {}
    self.require_all_exons_overlap = True
    self.output['overlap_length'] = 0
    #self.output['shared_exons'] = 0
    self.output['consecutive_exons'] = 0
    self.output['overlap_fractions'] = []
    #self.output['nonconsecutive_exons_matching'] = 0
    self.output['perfect_match'] = False
    self.output['full_match'] = False
    self.output['partial_match'] = False
    self.output['comparison_checked'] = False

  # Accept an overlap list of three floats, overlaps required for the first, internal and last exons
  def set_overlap_requirement(self,olist):
    self.overlap_requirement = olist    

  def set_require_all_exons_overlap(self,mybool):
    self.require_complete_overlap = mybool

  # Requires two genepred entries to compare
  def compare(self,eA,eB):
    self.output['comparison_checked'] = True
    range_list_A = [y[0] for y in sorted(eA.range_set.members, key=lambda x:x.get_payload())]
    first_exon_A = eA.get_first_exon_genomic_range()
    last_exon_A = eA.get_last_exon_genomic_range()
    range_list_B = [y[0] for y in sorted(eB.range_set.members, key=lambda x:x.get_payload())]
    first_exon_B = eB.get_first_exon_genomic_range()
    last_exon_B = eB.get_last_exon_genomic_range()
    self.output['overlap_length'] = 0
    best_consecutive = {}
    best_consecutive['exon_count'] = 0
    best_consecutive['overlap_size'] = 0
    best_consecutive['overlap_fractions'] = []
    current_consecutive = {}
    current_consecutive['exon_count'] = 0
    current_consecutive['overlap_size'] = 0
    current_consecutive['overlap_fractions'] = []

    ## Check for a totally perfect match
    #if len(range_list_A) == len(range_list_B):
    #  all_true = True
    #  for i in (0,len(range_list_A)):
    #    if not range_list_A[i].equals(range_list_B[i]):
    #      all_true = False
    #      break
    #  if all_true:
    #    print "found perfection"
    #    self.output['overlap_length'] = eA.length()
    #    self.output['consecutive_exons'] = len(range_list_A)
    #    self.output['overlap_fractions'] = [float(1) for x in range(0,len(range_list_A))]
    #    self.output['perfect_match'] = True
    #    self.output['full_match'] = True
    #    return

    ## Work on a nonperfect match
    indi = 0
    for rA in range_list_A:
      indi+=1
      any_count = 0
      match_size = 0
      match_count = 0 # depends on overlap
      match_fraction = 0
      indj = 0
      for rB in range_list_B:
        indj += 1

        # we only need to consider the same index exon for a 1 to 1 check
        if self.require_all_exons_overlap and indi != indj: continue
        if rA.equals(rB): # a perfect match is easy
          match_size += rA.length()
          match_fraction = 1
          match_count += 1
          any_count += 1
          continue
        reqover = self.overlap_requirement[1]
        if rA.equals(first_exon_A) or rB.equals(first_exon_B):
          reqover = self.overlap_requirement[0]
        elif rA.equals(last_exon_A) or rB.equals(last_exon_B):
          reqover = self.overlap_requirement[2]
        if not rA.overlaps(rB): continue
        osize = rA.overlap_size(rB)
        ofrac = float(osize)/rA.length()
        if rA.length() < rB.length(): ofrac = float(osize)/rB.length()
        any_count+=1 # count number of overlaping exons despite overlap
        if ofrac >= reqover:
          match_count += 1
          match_size += osize
          match_fraction = ofrac
        elif self.require_all_exons_overlap:
          # We are shortcutting if we have to have a full match to count anything
          return
      # Finished going through B
      if match_count > 1:
        break 

      #if its a mismatch then cut our 'best run saving'
      if match_count == 0:
        if current_consecutive['exon_count'] > best_consecutive['exon_count'] \
        or (current_consecutive['exon_count'] == best_consecutive['exon_count'] \
        and current_consecutive['overlap_size'] > best_consecutive['overlap_size']):
          best_consecutive['exon_count'] = current_consecutive['exon_count']
          best_consecutive['overlap_size'] = current_consecutive['overlap_size']
          best_consecutive['overlap_fractions'] = [x for x in current_consecutive['overlap_fractions']]
        current_consecutive['exon_count'] = 0
        current_consecutive['overlap_size'] = 0
        current_consecutive['overlap_fractions'] = []
      else:
        current_consecutive['exon_count'] += 1
        current_consecutive['overlap_size'] += match_size
        current_consecutive['overlap_fractions'].append(match_fraction)
    # Made it through A to no we can update our best one last time
    if current_consecutive['exon_count'] > best_consecutive['exon_count'] \
    or (current_consecutive['exon_count'] == best_consecutive['exon_count'] \
    and current_consecutive['overlap_size'] > best_consecutive['overlap_size']):
      best_consecutive['exon_count'] = current_consecutive['exon_count']
      best_consecutive['overlap_size'] = current_consecutive['overlap_size']
      best_consecutive['overlap_fractions'] = [x for x in current_consecutive['overlap_fractions']]

    # easy if theres zero consecutive matching exons
    if best_consecutive['exon_count'] == 0: return

    # If we're still here we must have some kind of match
    self.output['partial_match'] = True

    if best_consecutive['exon_count'] == len(range_list_A) \
    and best_consecutive['exon_count'] == len(range_list_B):
      self.output['full_match'] = True
    self.output['overlap_length'] = best_consecutive['overlap_size']
    self.output['consecutive_exons'] = best_consecutive['exon_count']
    self.output['overlap_fractions'] = [x for x in best_consecutive['overlap_fractions']]
    return

class GenePredFile:
  def __init__(self,filename):
    self.filename = filename
    self.gfr = GenericFileReader(filename)
    self.entries = []
    while True:
      line = self.gfr.readline()
      if not line: break
      if re.match('^#',line): continue
      gpe = GenePredEntry()
      gpe.line_to_entry(line)
      self.entries.append(gpe)
    return
  

# take an index-1 coordinate
# return true if it is present
def contains_coordinate(entry,coordinate):
  for i in range(0,len(entry['exonStarts'])):
    for j in range(entry['exonStarts'][i]+1,entry['exonEnds'][i]+1):
      if j == coordinate:
        return True
  return False

#pre: genepred entry, min gap size
#post: genepred entry with no gaps less than min gap size
def smooth_gaps(entry,gapsize):
  d = copy.deepcopy(entry)
  if d['exonCount'] == 1: return d
  prevEnd = d['exonEnds'][0]
  starts = []
  ends = []
  starts.append(d['exonStarts'][0])
  for i in range(1,d['exonCount']):
    gap = d['exonStarts'][i]-(prevEnd-1)+1
    if gap > gapsize:
      ends.append(prevEnd)
      starts.append(d['exonStarts'][i])
    prevEnd = d['exonEnds'][i]
  ends.append(prevEnd)
  d['exonStarts'] = starts
  d['exonEnds'] = ends
  d['exonCount'] = len(starts)
  if len(starts) != len(ends):
    sys.stderr.write("strange genepred error.\n")
    sys.exit()
  return d

#pre: genepred entry
#post: genepred line
#depreciated use entry_to_line
def genepred_entry_to_genepred_line(d):
  exonstarts = ",".join([str(x) for x in d['exonStarts']])+","
  exonends = ",".join([str(x) for x in d['exonEnds']])+","
  vals = [d['gene_name'], d['name'], d['chrom'], d['strand'], str(d['txStart']), str(d['txEnd']), str(d['cdsStart']), str(d['cdsEnd']), str(d['exonCount']),exonstarts,exonends]
  return "\t".join(vals)

#pre: genepred entry
#     0-based coordinate of the beginning of transcription
def left_trim_genepred(entry, coord):
  d = copy.deepcopy(entry)
  if coord <= d['txStart']: return d
  if coord+1 > d['txEnd']: 
    print "coordinate is out of bounds to left trim"
    sys.exit()
  exoncounts = 0
  starts = []
  ends = []  
  for i in range(0,d['exonCount']):
    if coord > d['exonEnds'][i]: 
      # this gets trimmed so skip ahead
      continue
    elif coord >= d['exonStarts'][i] and coord+1 <= d['exonEnds'][i]:
      #its in here, replace this
      exoncounts+=1
      starts.append(coord)
      ends.append(d['exonEnds'][i])
    else:
      # we have already trimmed
      exoncounts+=1
      starts.append(d['exonStarts'][i])
      ends.append(d['exonEnds'][i])
  d['exonStarts'] = starts
  d['exonEnds'] = ends
  d['exonCount'] = exoncounts
  d['txStart'] = d['exonStarts'][0]
  #see if we cut too much and need to extend back
  if d['txStart'] > coord:  d = left_extend_genepred(d,coord)
  return d

#pre: genepred entry
#     0-based coordinate of beginning of transcription
#post: new entry that lenthens the genepred to the coordinate prior to the start
def left_extend_genepred(entry,coord):
  d = copy.deepcopy(entry)
  if coord >= d['txStart']: return d
  d['txStart'] = coord
  d['exonStarts'][0] = coord
  return d

#pre: genepred entry
#     1-based coordinate of end of transcription
#post: new entry that lengthens the genepred to the coordinate beyond the end
def right_extend_genepred(entry,coord):
  d = copy.deepcopy(entry)
  if coord <= d['txEnd']: return d
  d['txEnd'] = coord
  d['exonEnds'][d['exonCount']-1] = coord
  return d

#pre: genepred entry
#     1-based coordinate of end of transcription
#post: new entry, that shortens the genepred to the coordinate
def right_trim_genepred(entry,coord):
  d = copy.deepcopy(entry)
  if coord >= d['txEnd']: return d
  if coord-1 < d['txStart']: 
    "coordinate is out of bounds for right trim"
    sys.exit()
  d['txEnd'] = coord
  exoncount = 0
  starts = []
  ends = []
  for i in range(0,d['exonCount']):
    exoncount += 1
    starts.append(d['exonStarts'][i])
    if d['exonEnds'][i] < coord:
      ends.append(d['exonEnds'][i])
      continue
      # breeze through
      # if we are still here we should replace
      # the end coordinate of this block, shortening it up
      # and then break out of here
    ends.append(coord)
    break
  d['exonStarts'] = starts
  d['exonEnds'] = ends
  d['exonCount'] = exoncount
  #see if we cut too much and need to extend back
  if d['txEnd'] < coord:  d = right_extend_genepred(d,coord)
  return d

# pre: chromosome
#      coordiante (1-indexed)
#      array of genepred entries for that chromosome
# post: return True if its the first base of an exon
def is_exon_start(chromosome, coordinate, genepred):
  for entry in genepred:
    if entry['chrom'] != chromosome:
      print "Error: you looking in the wrong chromosome of the genepred."
      sys.exit()
    if coordinate - 1 >= entry['txStart'] and coordinate <= entry['txEnd']:
      for exonStart in entry['exonStarts']:
        if exonStart == coordinate-1: return True
  return False

# pre: chromosome
#      coordiante (1-indexed)
#      array of genepred entries for that chromosome
# post: return True if its the last base of an exon
def is_exon_end(chromosome, coordinate, genepred):
  for entry in genepred:
    if entry['chrom'] != chromosome:
      print "Error: you looking in the wrong chromosome of the genepred."
      sys.exit()
    if coordinate - 1 >= entry['txStart'] and coordinate <= entry['txEnd']:
      for exonEnd in entry['exonEnds']:
        if exonEnd == coordinate: return True
  return False

def gene_annotate_by_coordinates(chromosome, start_index_0, end_index_1,annot_struct):
  if chromosome not in annot_struct:
    return
  genes = []
  for gene in annot_struct[chromosome]:
    i = 0
    for isoform in annot_struct[chromosome][gene]:
      i += 1
      geneset = [gene, isoform['strand'],isoform['txEnd']-isoform['txStart']]
      #encompassing
      if start_index_0 >= isoform['txStart'] and end_index_1 <= isoform['txEnd']:
        genes.append(geneset)
      elif start_index_0 <= isoform['txStart'] and end_index_1 >= isoform['txEnd']:
        genes.append(geneset)
      # left side
      elif start_index_0 > isoform['txStart'] and \
           start_index_0 <= isoform['txEnd']-1:
        genes.append(geneset)
      elif start_index_0 < isoform['txStart'] and \
           isoform['txStart'] <= end_index_1-1:
        genes.append(geneset)
      elif start_index_0 <= isoform['txEnd']-1 and \
           end_index_1 > isoform['txEnd']: 
        genes.append(geneset)
      elif isoform['txStart'] <= end_index_1-1 and \
           isoform['txEnd'] > end_index_1:
        genes.append(geneset)
  seen = set()
  norep = []
  for geneset in genes:
    if str(geneset) not in seen:
      norep.append(geneset)
    seen.add(str(geneset))

  #get non dash non Mir entries
  nondash = []
  dash = []
  for geneset in norep:
    p = re.compile('^MIR[\d]',re.IGNORECASE)
    if not re.search('-',geneset[0]) and not p.search(geneset[0]):
      nondash.append(geneset)
    else:
      dash.append(geneset)

  nondashsorted = sorted(nondash,key=lambda var: var[2], reverse=True)
  dashsorted = sorted(dash,key=lambda var: var[2], reverse=True)

  allsort = []
  for val in nondashsorted:
    allsort.append(val[0:2])
  for val in dashsorted:
    allsort.append(val[0:2])
  
  #remove nonredundant again
  seen = set()
  norep = []
  for geneset in allsort:
    if str(geneset) not in seen:
      norep.append(geneset)
    seen.add(str(geneset))

  return norep

# Pre: Genepred file name
# Post: A dictionary keyed by chromosome name with an array of all entries
def get_per_chromosome_array(genepred_filename):
  annot = {}
  i = 0
  with open(genepred_filename) as gpdfile:
    for line in gpdfile:
      if re.match('^#',line): continue
      d = line_to_entry(line)
      if d['chrom'] not in annot:
        annot[d['chrom']] = []
      annot[d['chrom']].append(d)
  return annot
        
# Pre: Genepred file name
# Post: A dictionary keyed by chromosome name with all genes on that chromosome there
#       with the txStart and txEnd for each isoform of each gene.
def get_gene_annotation_data_structure(genepred_filename):
  annot = {}
  i = 0
  with open(genepred_filename) as gpdfile:
    for line in gpdfile:
      if re.match('^#',line): continue
      d = line_to_entry(line)
      if d['chrom'] not in annot:
        annot[d['chrom']] = {}
      if d['gene_name'] not in annot[d['chrom']]:
        annot[d['chrom']][d['gene_name']] = []
      v = {}
      v['txStart'] = d['txStart']
      v['txEnd'] = d['txEnd']
      v['strand'] = d['strand']
      annot[d['chrom']][d['gene_name']].append(v)
  return annot
        
        
# Pre: Genepred file name
# Post: A dictionary keyed by transcript pointing to gene name
def get_transcript_to_gene_name_dictionary(genepred_filename):
  annot = {}
  with open(genepred_filename) as gpdfile:
    for line in gpdfile:
      if re.match('^#',line): continue
      d = line_to_entry(line)
      annot[d['name']] = d['gene_name']
  return annot
        
        

# pre: a genePred file
# post: a data structure for converting coordinates from transcripts described by a genepred
#       into reference coordinates 
#       dictionary keyed by sequence_name containing
#       'chrom' chromosome
#        or 'coordiantes' array of coordinates (zero indexed)
def get_directionless_gpd_conversion(genepred_filename):
  conv = {}
  with open(genepred_filename) as gpdfile:
    for line in gpdfile:
      if re.match('^#',line): continue
      d = line_to_entry(line)
      entry = {}
      entry['chrom'] = d['chrom']
      coords = []
      for i in range(0,d['exonCount']): 
        for j in range(d['exonStarts'][i],d['exonEnds'][i]): #account for that funky half open thing
          coords.append(j)
      entry['coordinates'] = coords
      conv[d['name']] = entry
  return conv

def entry_to_line(d):
  line = d['gene_name'] + "\t" + d['name'] + "\t" + d['chrom'] + "\t" \
       + d['strand'] + "\t" + str(d['txStart']) + "\t" + str(d['txEnd']) + "\t" \
       + str(d['cdsStart']) + "\t" + str(d['cdsEnd']) + "\t" +  str(d['exonCount']) +"\t" \
       + ','.join([str(x) for x in d['exonStarts']]) + ",\t" + ','.join([str(x) for x in d['exonEnds']]) + ','
  return line

# take a psl entry and a hash keyed on chromosome with the reference sequence
#      The latter is used to calculate the tSize
def entry_to_fake_psl_line(d,ref_hash):
  alen = 0
  qstarts = []
  for i in range(0,len(d['exonStarts'])):
    qstarts.append(alen)
    alen += d['exonEnds'][i]-d['exonStarts'][i]
  last = d['exonEnds'][len(d['exonEnds'])-1] #1-base
  first = d['exonStarts'][0] #0-base
  psl_line  = str(alen) + "\t"
  psl_line += "0\t0\t0\t0\t0\t0\t0\t" # misMatch through tBaseInsert
  psl_line += d['strand'] + "\t"
  psl_line += d['name'] + "\t"
  psl_line += str(last-first) + "\t"
  psl_line += "0\t"
  psl_line += str(last-first) + "\t"
  psl_line += d['chrom'] + "\t"
  psl_line += str(len(ref_hash[d['chrom']]))+"\t"
  psl_line += str(first) + "\t"
  psl_line += str(last) + "\t"
  psl_line += str(len(d['exonStarts'])) + "\t"
  psl_line += ','.join([str(d['exonEnds'][i]-d['exonStarts'][i]) for i in range(0,len(d['exonStarts']))])+','+"\t"
  psl_line += ','.join([str(x) for x in qstarts])+','+"\t"
  psl_line += ','.join([str(x) for x in d['exonStarts']])+','
  return psl_line

def line_to_entry(line):
  f = line.rstrip().split("\t")
  d = {}
  d['gene_name'] = f[0]
  d['name'] = f[1]
  d['chrom'] = f[2]
  d['strand'] = f[3]
  d['txStart'] = int(f[4])
  d['txEnd'] = int(f[5])
  d['cdsStart'] = int(f[6])
  d['cdsEnd'] = int(f[7])
  d['exonCount'] = int(f[8])
  exonstarts = [int(x) for x in f[9].rstrip(",").split(",")]
  d['exonStarts'] = exonstarts
  exonends = [int(x) for x in f[10].rstrip(",").split(",")]
  d['exonEnds'] = exonends
  return d

# pre: A genePred_file, a reference_fasta, an output fasta
# post: writes the output fasta file
#       negative strand transcripts are written in the positive orientation
#       this is done to make coordinate conversions easy for now
# modifies: file IO
def write_genepred_to_fasta_directionless(gpd_filename,ref_fasta,out_fasta):
  ofile = open(out_fasta,'w')
  ref = SequenceBasics.read_fasta_into_hash(ref_fasta)
  with open(gpd_filename) as f:
    for line in f:
      if re.match('^#',line): continue
      d = line_to_entry(line)
      if d['chrom'] in ref:
        seq = ''
        for i in range(0,d['exonCount']):
          seq = seq+ref[d['chrom']][d['exonStarts'][i]:d['exonEnds'][i]]
        ofile.write(">"+str(d['name'])+"\n"+seq+"\n")
  ofile.close()

# pre: A genePred_file, a reference_fasta, an output fasta
# post: writes the output fasta file
#        WARNING don't confuse this with the directionless version.
#        This version will reverse compelment outputs
# modifies: file IO
def write_genepred_to_fasta(gpd_filename,ref_fasta,out_fasta):
  ofile = open(out_fasta,'w')
  ref = SequenceBasics.read_fasta_into_hash(ref_fasta)
  with open(gpd_filename) as f:
    for line in f:
      if re.match('^#',line): continue
      d = line_to_entry(line)
      if d['chrom'] in ref:
        seq = ''
        for i in range(0,d['exonCount']):
          seq = seq+ref[d['chrom']][d['exonStarts'][i]:d['exonEnds'][i]]
        if d['strand'] == '-': seq = SequenceBasics.rc(seq)
        ofile.write(">"+str(d['name'])+"\n"+seq.upper()+"\n")
  ofile.close()

# pre: A genePred_file, an output genepred_file
# post: writes the output genepred file
#       where repeat names (field 1 index-0) have been renamed
# modifies: file IO
def write_uniquely_named_genepred(gpd_filename,gpd_out_filename):
  ofile = open(gpd_out_filename,'w')
  seennames = {}
  with open(gpd_filename) as f:
    for line in f:
      line = line.rstrip()
      if re.match('^#',line): 
        ofile.write(line+"\n")
        continue
      f = line.split("\t")
      if f[1] in seennames:
        seennames[f[1]] += 1
        f[1] = f[1] + "."+str(seennames[f[1]])
      else:
        seennames[f[1]] = 0
      ostring = "\t".join(f)
      ofile.write(ostring + "\n")
  ofile.close()

# Convert a bed file into genepred entries
# Uses bedtools merge so it throws away any depth information
# Pre: Min intron size, Max intron size, a sorted bed file
# Post: List of genepred entries
# Modifies:  Calls bedtools
def bed_to_genepred(min_intron_size,max_intron_size,bed_file):
  cmd = "bedtools merge -d "+str(min_intron_size)+" -i "+bed_file
  ps = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
  chrs = {}
  for line in ps.stdout:
    f = line.rstrip().split()
    chr = f[0]
    start = int(f[1])
    stop = int(f[2])
    if chr not in chrs: 
      chrs[chr] = {}
    chrs[chr][start]=stop
  ps.communicate()
  sorted_chrs = sorted(chrs.keys())
  gpds = []
  z = 0
  for chr in sorted_chrs:
    entry = []
    last_end = -1*max_intron_size*2
    starts = sorted(chrs[chr].keys())
    for start in starts:
      stop = chrs[chr][start]
      if start > last_end + max_intron_size:
        # output previous entry
        if len(entry) > 0:
          z+=1
          e = make_entry_from_starts_and_stops(str(z),str(z),chr,entry,'+')
          gpds.append(e)
        entry = []
      entry.append([start,stop])
      last_end = stop
    if len(entry) > 0:
      z+=1
      e = make_entry_from_starts_and_stops(str(z),str(z),chr,entry,'+')
      gpds.append(e)
  return gpds

def make_entry_from_starts_and_stops(name1,name2,chr,exons,strand):
  line = name1 + "\t" + name2 + "\t"
  line += chr + "\t"
  line += strand + "\t"
  line += str(exons[0][0]) + "\t" + str(exons[len(exons)-1][1]) + "\t"
  line += str(exons[0][0]) + "\t" + str(exons[len(exons)-1][1]) + "\t"
  line += str(len(exons)) + "\t"
  line += ",".join([str(x[0]) for x in exons])+"," + "\t"
  line += ",".join([str(x[1]) for x in exons])+","
  gpd = GenePredEntry()
  gpd.line_to_entry(line)
  return gpd  

# Pre: Take a filehandle a reading in GenePred entries that have been sorted by location
# Post: through the read_locus() method an array of genepred lines that are grouped by locus
#       these have not been grouped by strand or overlap, only by outer bounds
class GenePredLocusStream:
  def __init__(self,fhin):
    self.fh = fhin
    self.previous_line = self.fh.readline()
    self.finished = False
    self.range = None
    self.previous_range = None
    self.minimum_locus_gap = 0
    if not self.previous_line:
      self.finished = True
      return
    # initialize previous range
    gpd = GenePredEntry(self.previous_line)
    bed = RangeBasics.Bed(gpd.entry['chrom'],gpd.entry['txStart'],gpd.entry['txEnd'])
    self.previous_range = bed
    #print "initialize"
    return
  def set_minimum_locus_gap(self,ingap):
    self.minimum_locus_gap = ingap
  def read_locus(self):
    #print "open locus"
    if self.finished: return False
    buffer = []
    if self.previous_line:
      buffer.append(GenePredEntry(self.previous_line))
    while True:
      line = self.fh.readline()
      if not line:
        self.finished = True
        if len(buffer) > 0:
          return buffer
        return None
        # output buffer
      if self.different_locus(line): # We have finished one locus
        self.previous_line = line 
        #print 'outputing'
        return buffer
      buffer.append(GenePredEntry(line))
      self.previous_line = line

  def different_locus(self,gpd_line):
    gpd = GenePredEntry(gpd_line)
    bed = RangeBasics.Bed(gpd.entry['chrom'],gpd.entry['txStart'],gpd.entry['txEnd'])
    if not self.previous_range:
      self.previous_range = bed #update our range
      return True
    if bed.overlaps_with_padding(self.previous_range,self.minimum_locus_gap): # it overlaps with previous range
      #print 'new old range:'
      self.previous_range = self.previous_range.merge(bed)
      #print self.previous_range.get_range_string()
      return False
    self.previous_range = bed
    return True

# Pre: Take 2 filehandles reading in GenePreds that have been sorted by position
# Post: through the read_locus() method an array of genepred lines that are grouped by locus
#       these have not been grouped by strand or overlap, only by outer bounds
class GenePredDualLocusStream:
  def __init__(self,fhin1,fhin2):
    self.fh1 = fhin1
    self.fh2 = fhin2
    line1 = self.fh1.readline()
    line2 = self.fh2.readline()
    self.previous1 = None
    self.previous2 = None
    self.used1 = False
    self.used2 = False
    if line1:
      self.previous1 = GenePredEntry(line1)
    if line2:
      self.previous2 = GenePredEntry(line2)
    self.finished1 = False
    self.finished2 = False
    #self.range1 = None
    #self.range2 = None
    self.previous_range1 = None
    self.previous_range2 = None
    self.minimum_locus_gap = 0
    if not self.previous1:
      self.finished1 = True
    else:
      self.previous_range1 = self.previous1.get_bed()
      self.previous_range1.direction = None
    if not self.previous2:
      self.finished2 = True
    else:
      self.previous_range2 = self.previous2.get_bed()
      self.previous_range2.direction = None
    # initialize previous range
    return

  def set_minimum_locus_gap(self,ingap):
    self.minimum_locus_gap = ingap

  def read_locus(self):
    if self.finished1 and self.finished2: return False
    # Set our current locus to the lesser of the loaded in
    buffer = []
    buffer.append([])
    buffer.append([])
    buffer_range = None
    lesser = self.lesser_range()
    #print lesser
    if self.previous1 and (lesser == 1 or lesser == 0):
      buffer[0].append(self.previous1)
      buffer_range = self.previous1.get_bed()
      buffer_range.direction = None
    if self.previous2 and (lesser == 2 or lesser == 0):
      buffer[1].append(self.previous2)
      if not buffer_range: 
        buffer_range = self.previous2.get_bed()
        buffer_range.direction = None
      else: 
        b = self.previous2.get_bed()
        b.direction = None
        buffer_range = buffer_range.merge(b)
    done = False
    while True:
      if done or (self.finished1 and self.finished2):
        if len(buffer[0])==0 and len(buffer[1])==0: return None
        return buffer
      lesser = self.lesser_range() #1 or 2
      #print lesser
      overlaps = False
      if lesser == 1 or lesser == 0:
        # are about to use 1 so read in what will be the next one
        line1 = self.fh1.readline()
        if not line1:
          self.finished1 = True
          # output buffer
        else:
          g1 = GenePredEntry(line1)
          if not buffer_range: buffer_range = g1.get_bed()
          #check and see if we are overlapped
          if g1.get_bed().overlaps_with_padding(buffer_range,self.minimum_locus_gap):
            buffer[0].append(g1)
            b = self.previous_range1
            b.direction = None
            buffer_range = buffer_range.merge(b)
            overlaps = True
          #else: # we are done here
          #  done1 = True
          #  #print 'done1'
          self.previous1 = g1
          self.previous_range1 = g1.get_bed()
      if lesser == 2 or lesser == 0:
        line2 = self.fh2.readline()
        if not line2:
          self.finished2 = True
          # output buffer
        else:
          g2 = GenePredEntry(line2)
          if not buffer_range: buffer_range = g2.get_bed()
          if g2.get_bed().overlaps_with_padding(buffer_range,self.minimum_locus_gap):
            buffer[1].append(g2)
            b = self.previous_range2
            b.direction = None
            buffer_range = buffer_range.merge(b)     
            overlaps = True
          #else: # we are done here
          #  done2 = True
          #  #print 'done2'
          self.previous2 = g2
          self.previous_range2 = g2.get_bed()     
          self.previous_range2.direction = None
      if not overlaps: done = True
      #if done2 and lesser == 2:
      #  return buffer
      #if done1 and lesser == 1:
      #  return buffer
      #if self.different_locus(line1): # We have finished one locus
      #  self.previous_line1 = line1 
      #  #print 'outputing'
      #  return buffer1

  #def different_locus(self,gpd_line):
  #  gpd = GenePredEntry(gpd_line)
  #  bed = RangeBasics.Bed(gpd.entry['chrom'],gpd.entry['txStart'],gpd.entry['txEnd'])
  #  if not self.previous_range1:
  #    self.previous_range1 = bed #update our range
  #    return True
  #  if bed.overlaps_with_padding(self.previous_range1,self.minimum_locus_gap): # it overlaps with previous range
  #    #print 'new old range:'
  #    self.previous_range1 = self.previous_range1.merge(bed)
  #    #print self.previous_range.get_range_string()
  #    return False
  #  self.previous_range1 = bed
  #  return True

  def lesser_range(self):
    if self.finished1: return 2
    if self.finished2: return 1
    #print '--------------'
    #print self.previous_range1.get_range_string()
    #print self.previous_range2.get_range_string()
    c = self.previous_range1.cmp(self.previous_range2,self.minimum_locus_gap)
    if c == 0: return 0
    if c == -1: return 1
    else: return 2
