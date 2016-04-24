#!/usr/bin/python
import sys, argparse, StringIO, re, gzip
from multiprocessing import Pool, cpu_count, Queue, Lock
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
glock = Lock()
of = None
def main():
  global blocks
  parser = argparse.ArgumentParser(description="Take a bgzf compressed fastq file and make an index",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input_file',help="BGZF compressed fastq file")
  parser.add_argument('--threads',type=int,default=1,help="number of threads")
  args = parser.parse_args()
  if not is_bgzf(args.input_file):
    sys.stderr.write("ERROR: not a proper BGZF compressed file\n")
    sys.exit()
  sys.stderr.write("scanning block starts\n")
  bs = get_block_bounds(args.input_file)
  blocks[bs[0][0]] = [[bs[0][1],-1]]
  sys.stderr.write("scanning for new lines\n")
  z = 0
  scan_for_newlines(bs,args)
  # now need to clean out non start lines
  ncount = 0
  bnames = sorted(blocks.keys())
  for bname in bnames:
     v = []
     #print blocks[bname]
     for i in range(len(blocks[bname])):
       if ncount%4 == 0: v.append(blocks[bname][i])
       ncount += 1
     blocks[bname] = v
  sys.stderr.write("Traverse blocks and writing index\n")
  z = 0
  traverse_and_write(blocks,args)

def traverse_and_write(blocks,args):
  global of
  n = 0
  z = 0
  of = gzip.open(args.input_file+'.bgi','w')
  bnames = sorted(blocks.keys())
  chunk_size = 1000
  for block_set in [bnames[x:x+chunk_size] for x in range(0,len(bnames),chunk_size)]:
    ready_set = [[x,blocks[x][0][0],[y[1] for y in blocks[x]]] for x in block_set if len(blocks[x]) > 0]
    #block set is a block, 
    sys.stderr.write(str(n)+'/'+str(len(bnames))+"\r")
    n += chunk_size
    z += 1
    v = prepare_blocks(ready_set,args,z)
    write_results(v)

  sys.stderr.write("\n")
  of.close()

def write_results(v):
  global of
  for (i,line) in v:
    of.write(line)
  return

def prepare_blocks(blocks,args,z):
  results = []
  for block in blocks:
    #if len(blocks[block]) == 0: continue
    bstart = block[0]
    bend = block[1]
    starts = [x+1 for x in block[2]]
    #print starts
    with open(args.input_file,'rb') as inf:
      inf.seek(bstart)
      bytes = inf.read(bend-bstart)
      s = StringIO.StringIO(bytes)
      v = BGZF_reader(s)
      ubytes = v.read(70000)
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
          output = m.group(1)[1:]+"\t"+str(bstart)+"\t"+str(starts[i])+"\t"+str(len(m.group(1))+len(m.group(2))+len(m.group(3))+2)+"\t"+str(len(m.group(2)))+"\n"
          results.append([z,output])
    with open(args.input_file,'rb') as inf:
      v2 = BGZF_reader(inf,blockStart=bstart,innerStart=starts[-1]-1)
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
      output = m.group(1)[1:]+"\t"+str(cur)+"\t"+str(inn)+"\t"+str(len(buffer))+"\t"+str(len(m.group(2)))+"\n"
      results.append([z,output])
  return results

def scan_for_newlines(bs,args):
  z = 0
  # do them in chunks
  #for bounds in bs:
  chunk_size = 1000
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for bounds_blocks in [bs[x:x+chunk_size] for x in range(0,len(bs),chunk_size)]:
    if args.threads > 1:
      p.apply_async(get_nls,args=(bounds_blocks,args.input_file,z),callback=do_nls_output)
      #do_nls_output(v)
    else:
      v = get_nls(bounds_blocks,args.input_file,z)
      do_nls_output(v)
    sys.stderr.write(str(z)+'/'+str(len(bs))+"\r")
    z += chunk_size
  sys.stderr.write("\n")
  if args.threads > 1:
    p.close()
    p.join()

def get_nls(bounds_blocks,fname,i):
  breaks = []  
  for bounds in bounds_blocks:
    with open(fname,'rb') as inf:
      inf.seek(bounds[0])
      bytes = inf.read(bounds[1]-bounds[0])
      s = StringIO.StringIO(bytes)
    #v = BGZF_reader(inf,blockStart=bound[0],innerStart=0)
    v = BGZF_reader(s)
    ubytes = v.read(70000) # always less than 65K by definition
    p = re.compile('\n')
    nls = [m.start() for m in p.finditer(ubytes)]
    for j in range(len(nls)):
      breaks.append([bounds[0],bounds[1],nls[j],i])
    s.close()
  return breaks

def do_nls_output(results):
  global blocks
  #global ncount
  global glock
  glock.acquire()
  for e in results:
      #useval = False
      #if ncount%4 == 0: useval = True
      #ncount += 1
      #if not useval: continue
      if e[0] not in blocks: blocks[e[0]] = []
      blocks[e[0]].append([e[1],e[2]])
  glock.release()
    
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
