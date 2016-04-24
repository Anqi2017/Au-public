#!/usr/bin/python
import sys, argparse, StringIO, re, gzip
from multiprocessing import Pool, cpu_count, Queue
from Bio.Format.BGZF import is_bgzf, reader as BGZF_reader, get_block_bounds
from Bio.Format.Fastq import FastqEntry


# Create an index for bgzf zipped fastq files.  
# Pre: A fastq file that has been compressed by bgzf
# Post: the Pre file, with the exension .bgi added.
#       the index is gzipped
#       <name> <blockStart> <innerStart> <dataSize> <read length>
#       Be cautious that name could contain spaces and tabs

# global
blocks = {}
ncount = 1

def main():
  global blocks
  parser = argparse.ArgumentParser(description="Take a bgzf compressed fastq file and make an index",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input_file',help="BGZF compressed fastq file")
  parser.add_argument('--threads',type=int,default=1,help="number of threads")
  args = parser.parse_args()
  if not is_bgzf(args.input_file):
    sys.stderr.write("ERROR: not a proper BGZF compressed file\n")
    sys.exit()
  z = 0
  sys.stderr.write("scanning block starts\n")
  bs = get_block_bounds(args.input_file)
  blocks[bs[0][0]] = [[bs[0][1],-1]]
  sys.stderr.write("scanning for new lines\n")
  z = 0
  #if args.threads > 1:
  #  p = Pool(processes=args.threads)
  #results = []
  #for xs in [bs[j:j+args.threads*10] for j in range(0,len(bs),args.threads*10)]:
  for bounds in bs:
    #print xs
    #for bounds in xs:
    z += 1
    #if args.threads > 1:
    #  nls = p.apply_async(get_nls,args=(xs,args.input_file,z,))
    #else:
    #  nls = Queue()
    #  nls.put(get_nls(xs,args.input_file,z))
    v = get_nls(bounds,args.input_file,z)
    do_nls_output(v)
    #results.append(nls)
    sys.stderr.write(str(z)+'/'+str(len(bs))+"\r")
    #sys.exit()
  #if args.threads > 1:
  #  p.close()
  #  p.join()
  sys.stderr.write("\n")
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
      #print ubytes[starts[-2]:]
      for i in range(len(starts)-1):
        if starts[i] >= len(ubytes): #problem
          sys.stderr.write("Problem start\n")
          sys.exit()
        m = re.match('([^\n]+)\n([^\n]+)(\n[^\n]+\n[^\n]+)',ubytes[starts[i]:])
        if not m:
          sys.stderr.write("Problem overlap\n")
          sys.exit()
        else:
          if m.group(1)[0] != '@':
            sys.stderr.write("failed to parse last\n")
            sys.exit()  
          of.write(m.group(1)[1:]+"\t"+str(block)+"\t"+str(starts[i])+"\t"+str(len(m.group(1))+len(m.group(2))+len(m.group(3))+2)+"\t"+str(len(m.group(2)))+"\n")
    with open(args.input_file,'rb') as inf:
      v2 = BGZF_reader(inf,blockStart=block,innerStart=starts[-1]-1)
      spc = v2.read(1)
      if spc != "\n": 
        sys.stderr.write("expected newline\n")
        sys.exit()
      cur = v2.get_block_start()
      inn = v2.get_inner_start()
      buffer = ''
      for i in range(0,4):
        while True:
          c = v2.read(1)
          if len(c) == 0: break
          buffer += c
          if c == "\n": break
      if buffer == "":
        break
      m = re.match('([^\n]+)\n([^\n]+)',buffer)
      if not m:
        sys.stderr.write("failed to parse last\n"+buffer+"\n")
        sys.exit()
      if m.group(1)[0] != '@':
        sys.stderr.write("failed to parse last\n"+buffer+"\n")
        sys.exit()  
      of.write(m.group(1)[1:]+"\t"+str(cur)+"\t"+str(inn)+"\t"+str(len(buffer))+"\t"+str(len(m.group(2)))+"\n")

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

def get_nls(bounds,fname,i):
  with open(fname,'rb') as inf:
    inf.seek(bounds[0])
    bytes = inf.read(bounds[1]-bounds[0])
    s = StringIO.StringIO(bytes)
  #v = BGZF_reader(inf,blockStart=bound[0],innerStart=0)
  v = BGZF_reader(s)
  ubytes = v.read(70000) # always less than 65K by definition
  p = re.compile('\n')
  nls = [m.start() for m in p.finditer(ubytes)]
  breaks = []  
  for j in range(len(nls)):
    breaks.append([bounds[0],bounds[1],nls[j]])
  return breaks

def do_nls_output(results):
  global blocks
  global ncount
  #local = {}
  #for y in [x for x in results]:
  #  local[y[0]] = y[1]
  #for i in sorted(local):
  #  for e in local[i]:
  for e in results:
      #print e
      #print ncount
      useval = False
      if ncount%4 == 0: useval = True
      ncount += 1
      if not useval: continue
      if e[0] not in blocks: blocks[e[0]] = []
      blocks[e[0]].append([e[1],e[2]])
  #only every fourth newline is a start
  #breaks = [breaks[i] for i in range(0,len(breaks),4)]
  #sys.stderr.write("Reducing to new lines indicating starts\n")
  #blocks = {}
  #for i in range(0,len(breaks),4):
  #  if breaks[i][0] not in blocks: blocks[breaks[i][0]] = []
  #  blocks[breaks[i][0]].append([breaks[i][1],breaks[i][2]])
    
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
