import Bio.Structure
from Bio.Range import GenomicRange

# This whole format is a subclass of the Transcript subclass
class GPD(Bio.Structure.Transcript):
  def __init__(self,gpd_line):
    self._entry = self._line_to_entry(gpd_line)
    self._line = gpd_line.rstrip()
    self.exons = []
    self.junctions = []
    self._direction = self.value('strand')
    self._gene_name = self.value('gene_name')
    self._transcript_name = self.value('name')
    for i in range(0,self.value('exonCount')):
      ex = Bio.Structure.Exon(GenomicRange(self.value('chrom'),self.value('exonStarts')[i]+1,self.value('exonEnds')[i]))
      self.exons.append(ex)
    if self.value('exonCount') > 1:
      for i in range(0,self.value('exonCount')-1):
        l = GenomicRange(self.value('chrom'),self.value('exonEnds')[i],self.value('exonEnds')[i])
        r = GenomicRange(self.value('chrom'),self.value('exonStarts')[i+1]+1,self.value('exonStarts')[i+1]+1)
        junc = Bio.Structure.Junction(l,r)
        junc.set_exon_left(self.exons[i])
        junc.set_exon_right(self.exons[i+1])
        self.junctions.append(junc)


  def __str__(self):
    return self.get_gpd_line()  

  #output the original gpd line
  # Overrides Structure.Transcript
  def get_gpd_line(self):
    return self._line

  def value(self,key):
    return self._entry[key]

  def _line_to_entry(self,line):
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
