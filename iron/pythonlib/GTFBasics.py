# GTFBasics.py
#
# A module for holding functions for handling psl files
# 
# These include:
# line_to_entry() - read a line from a psl file into a dictionary
#           also includes a new array of adjusted coordinates 
#           in the same coordiante system as a positive strand
#           for easier comparison of query alignments
#

import re, sys

class GTFFile:
  def __init__(self,filename):
    self.genes = {}
    self.transcripts = {}
    with open(filename) as inf:
      for line in inf:
        if re.match('^#',line): continue
        f = line.rstrip().split("\t")
        if f[2] != 'exon': continue
        e = line_to_entry(line)
        if not e: continue
        if 'gene_id' not in e['attributes']:
          sys.stderr.write("WARNING no gene_id attribute found for "+line+"\n"+str(e)+"\n")
          continue
        if 'transcript_id' not in e['attributes']:
          sys.stderr.write("WARNING no gene_id attribute found for "+line+"\n")
          continue
        gene_id = e['attributes']['gene_id']
        transcript_id = e['attributes']['transcript_id']
        if gene_id not in self.genes:
          self.genes[gene_id] = []
        self.genes[gene_id].append(e)
        if transcript_id not in self.transcripts:
          self.transcripts[transcript_id] = []
        self.transcripts[transcript_id].append(e)
    return
  def write_genepred(self,filehandle):
    for transcript in self.transcripts:
      #filehandle.write(str(transcript)+"\n")
      tlist = self.transcripts[transcript]
      elist = {}
      gene = '.'
      strand = '.'
      chrom = '.'
      for t in tlist:
        if not t['gff'][2].lower() == 'exon':
          continue
        start = int(t['gff'][3])-1
        end = int(t['gff'][4])
        if start > end: 
          sys.stderr.write("ERROR start bigger than end\n")
          sys.exit()
        elist[start] = end
        gene = t['attributes']['gene_id']
        strand = t['gff'][6]
        chrom = t['gff'][0]
      #print strand
      #print transcript
      starts = sorted(elist,key=lambda k: elist[k])
      first = starts[0]
      last = elist[starts[len(starts)-1]]
      ostring = gene + "\t" + transcript + "\t" + chrom + "\t" + strand + "\t" + str(first) + "\t" \
              + str(last) + "\t" + str(first) + "\t" + str(last) + "\t" \
              + str(len(starts)) + "\t" \
              + ",".join([str(x) for x in starts]) + ",\t" \
              + ",".join([str(elist[x]) for x in starts])+","
      filehandle.write(ostring+"\n")

class GTF:
  def __init__(self,line=None):
    self.entry = None
    if line: self.entry = line_to_entry(line)

def line_to_entry(line):
  f = line.rstrip().split("\t")
  gff_fields = f[:8]
  preattributes = re.split('\s*;\s*',f[8])
  attributes = {}
  for attribute in preattributes:
    m = re.search('(\S+)\s*["\']([^\'"]+)["\']',attribute)
    if m:  
      attributes[m.group(1)] = m.group(2)
  if len(attributes.keys()) > 0:
    entry = {}
    entry['gff'] = gff_fields
    entry['attributes'] = attributes
    return entry
  return None
