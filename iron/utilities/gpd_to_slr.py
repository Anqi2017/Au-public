#!/usr/bin/python
import argparse, sys, os, gzip, random
from shutil import rmtree
from multiprocessing import cpu_count, Lock, Pool
from tempfile import mkdtemp, gettempdir
from uuid import uuid4
from subprocess import PIPE, Popen

from Bio.Format.Fasta import FastaData
from Bio.Format.GPD import GPD
from Bio.Format.PSL import PSL
from Bio.Sequence import Seq

glock = Lock()
#manager = Manager()
of = None

def main():
  #do our inputs
  args = do_inputs()
  global of
  of = sys.stdout
  if args.output:
    if args.output[-4:] == '.bam':
      cmd = 'samtools view -Sb - -o '+args.output
      p = Popen(cmd.split(),stdin=PIPE)
      of = p.stdin
    else:
      sys.stderr.write("ERROR: stdout and .bam are the only valid output formats\n")
      sys.exit()
  inf = sys.stdin
  if args.input != '-':
    if args.input[-3:] == '.gz':
      inf = gzip.open(args.input)
    else: inf = open(args.input)
  sys.stderr.write("reading reference genome\n")
  ref = FastaData(open(args.reference).read())
  #shared = manager.dict()
  shared = {}
  for chr in sorted(ref.keys()): 
    sys.stderr.write("reading "+chr+"\n")
    shared[chr] = ref[chr].upper()
    ref.remove(chr)
  sys.stderr.write("finished reading shared memory reference\n")
  sys.stderr.write("Now make the header\n")
  of.write("@HD\tVN:1.0\tSO:unknown\n")
  of.write("@PG\tID:SLR\n")
  for chr in sorted(shared.keys()):
    of.write("@SQ\tSN:"+chr+"\tLN:"+str(len(shared[chr]))+"\n")

  if args.threads > 1:
    poo = Pool(processes=args.threads)

  buffer = []
  max_buffer = 1
  z = 0
  for line in inf:
    z += 1
    if z%1000==0: sys.stderr.write(str(z)+"   \r")
    buffer.append(line)
    if len(buffer) >= max_buffer:
      if args.threads == 1:
        results = do_buffer(buffer,shared,args)
        do_out(results)
      else:
        poo.apply_async(do_buffer,args=(buffer[:],shared,args,),callback=do_out)
      buffer = []
  if len(buffer) > 0:
    if args.threads ==1:
      results = do_buffer(buffer,shared,args)
      do_out(results)
    else:
      poo.apply_async(do_buffer,args=(buffer[:],shared,args,),callback=do_out)

  if args.threads > 1:
    poo.close()
    poo.join()

  sys.stderr.write("\n")
  if args.output:
    p.communicate()
  else: of.close()

  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_out(results):
  global glock
  glock.acquire()
  global of
  for result in results:
    of.write(result+"\n")
  glock.release()
  return

def do_buffer(gpd_lines,fasta,args):
  results = []
  for gpd_line in gpd_lines:
    gpd = GPD(gpd_line)
    l = gpd.get_length()
    if l < args.length: continue
    num = int(float(l)/float(args.length))
    rem = l % args.length
    #print 'rem : '+str(rem)
    extra = 0
    offset = 0
    #if space > 1: # we have room to make multiple passes
    #  #print '---'
    #  #print 'length: '+str(l)
    #  #print 'strand: '+gpd.get_strand()
    #  if random.random() < 0.5: extra = rem
    #  offset = int(float(args.length)/float(args.coverage))
    #else:
    #  offset = int(float(rem)/float(args.coverage)) 

    if args.short_reads:
      offset = 0
      if random.random() < 0.5: offset = rem
      gsub = gpd.subset(offset,args.length+offset)
      #print gsub.get_gpd_line()
      val = get_sam(gsub,fasta)
      results.append(val)
      #continue
    else:# not short reads
      for i in range(0,args.coverage):
        init = 0
        if num == 0 and rem > 0:
          init = random.choice(range(0,rem))
        elif num > 0:
          init = random.choice(range(0,args.length))
        #start = (i*offset+extra) % args.length
        #while start+args.length <= l:
        for j in range(init,l,args.length):
          if j + args.length > l: break
          #print str(start)+" "+str(start+args.length)
          gsub = gpd.subset(j,j+args.length)
          val = get_sam(gsub,fasta)
          results.append(val)
          #print gsub.get_sequence(fasta)
          #start += args.length
          #print gsub.get_strand()
    #print space
    #print rem
    #print gpd
  return results

def get_sam(gsub,fasta):
  gsub.set_transcript_name(str(uuid4()))
  psl = PSL(gsub.get_fake_psl_line(fasta))
  psl.set_reference = fasta
  seq = gsub.get_sequence(fasta)
  seq = Seq(seq).rc().seq
  psl.set_query_sequence(seq)
  sam = psl.get_SAM()
  return sam.get_line()

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('reference',help="Reference FASTA")
  parser.add_argument('-l','--length',type=int,default=100,help="length of output alignment")
  parser.add_argument('-c','--coverage',type=int,default=3,help="number of times to cover")
  parser.add_argument('--short_reads',action='store_true',help="instead of coverage, just sample l bases from left or right randomly once per read")  
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=1,help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  args = parser.parse_args()
  # Temporary working directory step 2 of 3 - Creation
  setup_tempdir(args)
  return args

def setup_tempdir(args):
  if args.specific_tempdir:
    if not os.path.exists(args.specific_tempdir):
      os.makedirs(args.specific_tempdir.rstrip('/'))
    args.tempdir = args.specific_tempdir.rstrip('/')
    if not os.path.exists(args.specific_tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  else:
    args.tempdir = mkdtemp(prefix="weirathe.",dir=args.tempdir.rstrip('/'))
    if not os.path.exists(args.tempdir.rstrip('/')):
      sys.stderr.write("ERROR: Problem creating temporary directory\n")
      sys.exit()
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

if __name__=="__main__":
  main()
