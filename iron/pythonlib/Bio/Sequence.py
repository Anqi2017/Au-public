import re, sys, base64, zlib
from string import maketrans
#from multiprocessing.sharedctypes import RawArray

class FastqHandleReader:
  def __init__(self,in_handle):
    self.fh = in_handle
  def close(self):
    self.fh.close()
  def read_entry(self):
    t = {}
    t['name'] = ''
    t['seq'] = ''
    t['original'] = ''
    t['qual'] = ''
    line1 = self.fh.readline()
    if not line1: return None
    line2 = self.fh.readline()
    if not line2: return None
    line3 = self.fh.readline()
    if not line3: return None
    line4 = self.fh.readline()
    if not line4: return None
    # end of the line, then finish it
    m1 = re.match('^@(.*)$',line1.rstrip())
    if not m1:
      sys.stderr.write('ERROR: '+str(line1)+"\n")
      sys.exit()
    t['name'] = m1.group(1)
    t['seq'] = line2.rstrip()
    t['qual'] = line4.rstrip()
    t['original'] = line1+line2+line3+line4
    return t
# an upgrade to the old sequence_basics set

    
def rc(seq):
  complement = maketrans('ACTGUNXactgunx','TGACANXtgacanx')
  return seq.translate(complement)[::-1]

# pre: A coordiante array
# post: a more bed format like string
def collapse_coordinate_array(readcoords):
  sc = readcoords[0]
  lc = readcoords[0] 

  # readcoodinates should be size of read length unless a sequnece in the read didn't map (was an insertion)
  cstring = ''
  for i in range(1,len(readcoords)):
    if readcoords[i] > int(lc)+1:
      cstring = cstring + ',' + str(sc) + "-" + str(lc)
      sc = readcoords[i]
    lc = readcoords[i]
  cstring = cstring + ',' +  str(sc) + "-" + str(lc)
  cstring = cstring.lstrip(',')
  return cstring

def encode_name(conversion_string):
  compressed_string = zlib.compress(conversion_string,9)
  enc_string = base64.b32encode(compressed_string)
  return 'SZ_'+enc_string.rstrip('=')

def decode_name(safename):
  frag = safename.lstrip('SZ_')
  padding = (8-(len(frag) % 8)) % 8
  c = base64.b32decode(frag+'='*padding)
  return  zlib.decompress(c)

