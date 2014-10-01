import genepred_basics, sequence_basics
import re

# pre: sam file name, genepred name 
# post: list of coordinates in the reference format
#       these coordinates are zero indexed for both the start and end coordinate
#       note thats different than gpd or psl or even wig or bed
#       <read name> <genepred entry name> <chromosome:coord1-coord2,coord3-coord4,...>
def convert_directionless_gpd_alignment_to_reference(sam_filename,genepred_filename,out_map):
  conv = genepred_basics.get_directionless_gpd_conversion(genepred_filename)
  ofile = open(out_map,'w')
  with open(sam_filename) as samfile:
    for line in samfile:
      line = line.rstrip()
      if re.match('^@[A-Z][A-Z]\s',line): continue #skip header
      d = sam_line_to_dictionary(line)
      if d['rname'] == '*': continue #skip unmapped
      startposition = d['pos']-1
      readcoord = []
      z = 0
      for entry in d['cigar_array']:
        if re.match('[MISX=]',entry['op']):  # all the entries that map to the read
          for i in range(0,entry['val']):
            if re.match('[M=X]',entry['op']): #all the entries that match the reference alignment
              readcoord.append(conv[d['rname']]['coordinates'][startposition+z])
              z+=1
            # lets ignore insertions for now
            #else:
            #  readcoord.append('*')
        if re.match('[DNH]',entry['op']):
          z+= entry['val']      
      abbrev = conv[d['rname']]['chrom']+':'+sequence_basics.collapse_coordinate_array(readcoord)
      ofile.write(d['qname'] + "\t" + d['rname'] + "\t" + abbrev + "\n")
  ofile.close()

# pre:       A line from a sam file
# post:      a string with the coordiantes of the alignment

def get_coordinates(sam_line):
  f = sam_line.rstrip().split("\t")
  name = f[0]
  coordinate = ''
  if f[2] == '*':
    return [name, coordinate]
  coordinate = f[2]+':'+str(f[3])+':'+f[5]
  return [name,coordinate]


#pre: a flag from a sam file, in integer format
#     a bit to convert, given as a hex number ie 0x10
#post: returns true if the flag is set on
def is_header(line):
  if re.match('^@',line):
    f = line.rstrip().split("\t")
    if(len(f) > 9):
      return False
    return True
  return False


#pre: a flag from a sam file, in integer format
#     a bit to convert, given as a hex number ie 0x10
#post: returns true if the flag is set on
def check_flag(flag,inbit):
  if flag & inbit:
    return True
  return False

#pre: a line from a sam file that is not a header entry
#post: a dictionary with entries named like the manual
def sam_line_to_dictionary(line):
  f = line.rstrip().split("\t")
  d = {}
  d['qname'] = f[0]
  d['flag'] = int(f[1])
  d['rname'] = f[2]
  d['pos'] = int(f[3])
  d['mapq'] = f[4]
  cigar = parse_cigar(f[5])
  d['cigar_array'] = cigar
  d['cigar'] = f[5]
  d['rnext'] = f[6]
  d['pnext'] = int(f[7])
  d['tlen'] = int(f[8])
  d['seq'] = f[9]
  d['qual'] = f[10]
  return d

# pre: CIGAR string
# post: an array of cigar string entries
def parse_cigar(cigar):
  v = re.findall('([0-9]+[A-Z])',cigar)
  vals = []
  for val in v:
    m = re.match('(\d+)([A-Z])',val)
    d = {}
    d['op'] = m.group(2)
    d['val'] = int(m.group(1))
    vals.append(d)
  return vals
