import re, sys, copy
import sequence_basics

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
    print "strange genepred error."
    sys.exit()
  return d

#pre: genepred entry
#post: genepred line
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
      d = genepred_line_to_dictionary(line)
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
      d = genepred_line_to_dictionary(line)
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
      d = genepred_line_to_dictionary(line)
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
      d = genepred_line_to_dictionary(line)
      entry = {}
      entry['chrom'] = d['chrom']
      coords = []
      for i in range(0,d['exonCount']): 
        for j in range(d['exonStarts'][i],d['exonEnds'][i]): #account for that funky half open thing
          coords.append(j)
      entry['coordinates'] = coords
      conv[d['name']] = entry
  return conv

def genepred_line_to_dictionary(line):
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
  ref = sequence_basics.read_fasta_into_hash(ref_fasta)
  with open(gpd_filename) as f:
    for line in f:
      if re.match('^#',line): continue
      d = genepred_line_to_dictionary(line)
      if d['chrom'] in ref:
        seq = ''
        for i in range(0,d['exonCount']):
          seq = seq+ref[d['chrom']][d['exonStarts'][i]:d['exonEnds'][i]]
        ofile.write(">"+str(d['name'])+"\n"+seq.upper()+"\n")
  ofile.close()

# pre: A genePred_file, a reference_fasta, an output fasta
# post: writes the output fasta file
#        WARNING don't confuse this with the directionless version.
#        This version will reverse compelment outputs
# modifies: file IO
def write_genepred_to_fasta(gpd_filename,ref_fasta,out_fasta):
  ofile = open(out_fasta,'w')
  ref = sequence_basics.read_fasta_into_hash(ref_fasta)
  with open(gpd_filename) as f:
    for line in f:
      if re.match('^#',line): continue
      d = genepred_line_to_dictionary(line)
      if d['chrom'] in ref:
        seq = ''
        for i in range(0,d['exonCount']):
          seq = seq+ref[d['chrom']][d['exonStarts'][i]:d['exonEnds'][i]]
        if d['strand'] == '-': seq = sequence_basics.rc(seq)
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

