#!/usr/bin/python
import sys, argparse, StringIO, re, gzip
from Bio.Format.BGZF import is_bgzf, reader as BGZF_reader, get_block_bounds
from Bio.Format.Fastq import FastqEntry

def main():
  parser = argparse.ArgumentParser(description="Take a bgzf compressed fastq file and make an index",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input_file',help="BGZF compressed fastq file")
  args = parser.parse_args()
  if not is_bgzf(args.input_file):
    sys.stderr.write("ERROR: not a proper BGZF compressed file\n")
    sys.exit()
  z = 0
  sys.stderr.write("scanning block starts\n")
  bs = get_block_bounds(args.input_file)
  breaks = [[bs[0][0],bs[0][1],-1]]
  sys.stderr.write("scanning for new lines\n")
  for i in range(0,len(bs)):
    nls = get_nls(bs[i],args.input_file)
    for j in range(len(nls)):
      breaks.append([bs[i][0],bs[i][1],nls[j]])
  #only every fourth newline is a start
  #breaks = [breaks[i] for i in range(0,len(breaks),4)]
  blocks = {}
  for i in range(0,len(breaks),4):
    if breaks[i][0] not in blocks: blocks[breaks[i][0]] = []
    blocks[breaks[i][0]].append([breaks[i][1],breaks[i][2]])
  sys.stderr.write("Traverse blocks and writing index\n")
  of = gzip.open(args.input_file+'.bgi','w')
  z = 0
  for block in sorted(blocks):
    z+=1
    sys.stderr.write(str(z)+'/'+str(len(blocks))+"\r")
    if len(blocks[block]) == 0: continue
    bend = blocks[block][0][0]
    starts = [x[1]+1 for x in blocks[block]]
    with open(args.input_file,'rb') as inf:
      inf.seek(block)
      bytes = inf.read(bend-block)
      s = StringIO.StringIO(bytes)
      v = BGZF_reader(s)
      ubytes = v.read(70000)
      # now we can find all the new starts
      # do all but the last
      for i in range(len(starts)-1):
        if starts[i] >= len(ubytes): #problem
          sys.stderr.write("Problem start\n")
          sys.exit()
        m = re.match('([^\n]+)\n',ubytes[starts[i]:])
        if not m:
          sys.stderr.write("Problem overlap\n")
          sys.exit()
        else:
          of.write(m.group(1)+"\t"+str(block)+"\t"+str(starts[i])+"\n")
    with open(args.input_file,'rb') as inf:
      v2 = BGZF_reader(inf,blockStart=block,innerStart=starts[-1]-1)
      spc = v2.read(1)
      if spc != "\n": 
        sys.stderr.write("expected newline\n")
        sys.exit()
      cur = v2.get_block_start()
      inn = v2.get_inner_start()
      buffer = ''
      while True:
        c = v2.read(1)
        if len(c) == 0: break
        if c == "\n": break
        buffer += c
      of.write(buffer+"\t"+str(cur)+"\t"+str(inn)+"\n")

  sys.stderr.write("\n")
  sys.exit()
  buffer = ''     
  with open(args.input_file) as inf:
    #inf.seek(bs[i])
    reader = BGZF_reader(inf)
    while True:
        cur = reader.get_block_start()
        inn = reader.get_inner_start()
        fq = readfastq(reader)
        z += 1
        if not fq: break
        if z%1000 == 0: sys.stderr.write("Indexed "+str(z)+" reads\r")
        of.write(fq['name']+"\t"+str(cur)+"\t"+str(inn)+"\n")
    inf.close()
    sys.stderr.write("\n")
  of.close()

def get_nls(bounds,fname):
  with open(fname,'rb') as inf:
    inf.seek(bounds[0])
    bytes = inf.read(bounds[1]-bounds[0])
    s = StringIO.StringIO(bytes)
  #v = BGZF_reader(inf,blockStart=bound[0],innerStart=0)
  v = BGZF_reader(s)
  ubytes = v.read(70000) # always less than 65K by definition
  p = re.compile('\n')
  return [m.start() for m in p.finditer(ubytes)]
    
def readfastq(reader):
  buffer = ''
  for i in range(0,4):
    v = ''
    while v!="\n":
      v = reader.read(1)
      if len(v) == 0: return False
      buffer += v
  if len(buffer) == 0: return False
  return FastqEntry(buffer.rstrip().split("\n"))

if __name__=="__main__":
  main()
