#!/usr/bin/python
import argparse, sys, os
from shutil import rmtree, move
from multiprocessing import cpu_count, Pool, Lock, Queue
from tempfile import mkdtemp, gettempdir
import math
import re


##################################
# A python parallel appraoch to
# mergesort

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Merge sort. Low memory multi-threaded.  Not as fast as unix sort.")
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT if not set")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
  group.add_argument('--memory','-m',action='store_true',help="Do sort in memory")
  parser.add_argument('--buffer_size',default=100000,type=int,help="INT Number of lines to sort at at time")
  parser.add_argument('--fields','-f',help="Search fields '1,2n,3ni' would do field 1 first, then field 2 numerically,then field 3 numerically but inverted")
  parser.add_argument('--maxbytes',type=int,default=100000000,help="Max temporary file to pull into memory\n")
  args = parser.parse_args()
  # Setup inputs 
  if args.input == '-':
    args.input = sys.stdin
  else:
    args.input = open(args.input)
  # Temporary working directory step 2 of 3 - Creation
  if not args.memory:
    setup_tempdir(args)
  return args

def main():
  #do our inputs
  args = do_inputs()
  #1. Stream through the initial sorts
  buffer = []
  cnt = 0
  results = [] #only filled if we are using args.memory
  if args.threads > 1:
    p = Pool(processes=args.threads)
  for line in args.input:
    buffer.append(line)
    if len(buffer) >= args.buffer_size:
      cnt += 1
      if args.threads > 1:
        r = p.apply_async(process_buffer,args=(buffer,cnt,args))
        results.append(r)
      else:
        r = process_buffer(buffer,cnt,args) #consider making copy of buffer here with slice but i think it gets a new copy anyways on multiprocessing
        q = Queue()
        q.put(r)
        results.append(q)
      buffer = []
  if len(buffer) > 0: 
    cnt += 1
    if args.threads > 1:
      r = p.apply_async(process_buffer,args=(buffer,cnt,args))
      results.append(r)
    else:
      r = process_buffer(buffer,cnt,args)
      q = Queue()
      q.put(r)
      results.append(q)
  if args.threads > 1:
    p.close()
    p.join()
  #2. merge the files from the bottom up
  if not args.memory:
    while cnt != 1:
      cnt = bottom_up(args,cnt)
    #3. Do output
    # Setup outputs
    if args.output:
      args.output = open(args.output,'w')
    else:
      args.output = sys.stdout
    with open(args.tempdir+'/l.1') as inf:
      for line in inf:
        args.output.write(line)
    args.output.close()
  else: # do the memory way
    while len(results) > 1:
      results = bottom_up_mem(results,args)    
    #for r in results:
    #v = list(results[0].get())
    #sys.exit()
    if args.output:
      args.output = open(args.output,'w')
    else:
      args.output = sys.stdout
    for vals in results:
      for val in vals.get():
        args.output.write(val)
    args.output.close()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir and not args.memory:
    rmtree(args.tempdir)

def bottom_up_mem(results,args):
  rlen = len(results)
  newresults = []
  if args.threads > 1:
    p = Pool(processes=args.threads)
  while(len(results) > 0):
    r1t = results.pop(0)
    r1 = r1t.get()
    r2 = None
    if len(results) > 0:
      r2t = results.pop(0)
      r2 = r2t.get()
    if args.threads > 1:
      r = p.apply_async(merge_mem,args=(r1,r2,args))
      newresults.append(r)
    else:
      r = merge_mem(r1,r2,args)
      q = Queue()
      q.put(r)
      newresults.append(q)
  #now lets move the names back to 'l' type to make recursion easy
  if args.threads > 1:
    p.close()
    p.join()
  return newresults

def do_compare(line1,line2,args):
  if not args.fields:
    if line1 < line2:
      return True
    return False
  fields = args.fields.split(',')
  lf1 = line1.rstrip().split("\t")
  lf2 = line2.rstrip().split("\t")
  for f in fields:
    myTrue = True
    myFalse = False
    #see if we are inverting
    if re.search('i',f):
      myTrue = False
      myFalse = True
    m = re.search('(\d+)',f)
    if not m:
      sys.stderr.write("ERROR: must specify a field index (base-1) to sort on sort on with fields option\n")
      sys.exit()
    i = int(m.group(1))-1
    isN = False
    if re.search('n',f): isN = True
    if isN:
      if float(lf1[i]) < float(lf2[i]):
        return myTrue
      elif float(lf1[i]) > float(lf2[i]):
        return myFalse
    else:
      if lf1[i] < lf2[i]:
        return myTrue
      elif lf1[i] > lf2[i]:
        return myFalse
  if line1 < line2:  #default to string sort
    return True
  return False

def merge_mem(r1,r2,args):
  # case where f2 is not there is first
  if not r2:
    return r1
  # case where we merge by file
  r1ind = [0]
  r2ind = [0]
  r1len = len(r1)
  r2len = len(r2)
  used1 = True
  used2 = True
  inf1 = None
  inf2 = None
  line1 = get_line(inf1,r1,r1ind,r1len)
  line2 = get_line(inf2,r2,r2ind,r2len)
  rout = []
  while True:
    if not line1 and not line2:
      break #at both EOFs
    if line1 and not line2:
      rout.append(line1)
      # finish 1
      while True:
        line1 = get_line(inf1,r1,r1ind,r1len)
        if not line1: break
        rout.append(line1)
      break        
    elif line2 and not line1:
      rout.append(line2)
      # finish 2
      while True:
        line2 = get_line(inf2,r2,r2ind,r2len)
        if not line2: break
        rout.append(line2)
      break
    elif do_compare(line1,line2,args):
      rout.append(line1)
      line1 = get_line(inf1,r1,r1ind,r1len)
    else:
      rout.append(line2)
      line2 = get_line(inf2,r2,r2ind,r2len)
  # Finished merge now clean up
  return rout

def merge_files(f1,f2,cnt,args):
  # case where f2 is not there is first
  if not f2:
    move(f1,args.tempdir+'/m.'+str(cnt))
    return
  # case where we merge by file
  of = open(args.tempdir+'/m.'+str(cnt),'w')
  f1size = os.path.getsize(f1)
  f2size = os.path.getsize(f2)
  useFiles = True
  inf1 = None
  inf2 = None
  f1lines =None
  f2lines = None
  f1ind = [0]
  f2ind = [0]
  f1len = None
  f2len = None
  if f1size <= args.maxbytes and f2size <= args.maxbytes:
    useFiles = False
  if not useFiles: #get arrays ready if the data is small
    f1lines = []
    f2lines = []
    # do them all in memory
    with open(f1) as inf:
      for line in inf: f1lines.append(line)
    f1len = len(f1lines)
    with open(f2) as inf:
      for line in inf: f2lines.append(line)
    f2len = len(f2lines)
  else:
    inf1 = open(f1)
    inf2 = open(f2)
  used1 = True
  used2 = True
  line1 = get_line(inf1,f1lines,f1ind,f1len)
  line2 = get_line(inf2,f2lines,f2ind,f2len)
  #line1 = inf1.readline()
  #line2 = inf2.readline()
  while True:
    if not line1 and not line2:
      break #at both EOFs
    if line1 and not line2:
      of.write(line1)
      # finish 1
      while True:
        line1 = get_line(inf1,f1lines,f1ind,f1len)
        if not line1: break
        of.write(line1)
      break        
    elif line2 and not line1:
      of.write(line2)
      # finish 2
      while True:
        line2 = get_line(inf2,f2lines,f2ind,f2len)
        if not line2: break
        of.write(line2)
      break
    elif do_compare(line1,line2,args):
      of.write(line1)
      line1 = get_line(inf1,f1lines,f1ind,f1len)
    else:
      of.write(line2)
      line2 = get_line(inf2,f2lines,f2ind,f2len)
  # Finished merge now clean up
  if useFiles:
    inf1.close()
    inf2.close()
  of.close()
  os.remove(f1)
  os.remove(f2)

def get_line(fh,farr,find,flen):
  if not fh and not farr:
    return None
  if fh:
    return fh.readline()
  if find[0] < flen:
    curr = find[0]
    find[0]+=1
    return farr[curr]
  return None

def bottom_up(args,cnt):
  if args.threads > 1:
    p = Pool(processes=args.threads)
  files = [args.tempdir+'/l.'+str(x+1) for x in range(0,cnt)]
  newcount = 0
  while(len(files) > 0):
    f1 = files.pop(0)
    f2 = None
    if len(files) > 0:
      f2 = files.pop(0)
    newcount += 1
    if args.threads > 1:
      p.apply_async(merge_files,args=(f1,f2,newcount,args))
    else:
      merge_files(f1,f2,newcount,args)
  #now lets move the names back to 'l' type to make recursion easy
  if args.threads > 1:
    p.close()
    p.join()
  for i in range(1,newcount+1):
    move(args.tempdir+'/m.'+str(i),args.tempdir+'/l.'+str(i))
  return newcount

def process_buffer(buffer,cnt,args):
  sorted_buffer = merge_sort(buffer,args)
  if args.memory:
    # we just need the buffer if we're staying in memory lane
    return sorted_buffer
  of = open(args.tempdir+"/l."+str(cnt),'w')
  for line in sorted_buffer:
    of.write(line)
  of.close()
  return None

  #sys.stderr.write(str(len(sorted_buffer))+"\n")
  return [sorted_buffer,cnt,args]

def write_temp(s):
  return 

def merge_sort(a,args):
  length_a = len(a)
  if length_a <= 1: return a
  m = int(math.floor(length_a/2))
  a_left = a[0:m]
  a_right = a[m:]
  a_left = merge_sort(a_left,args)
  a_right = merge_sort(a_right,args)
  return merge(a_left,a_right,args)

def merge(left,right,args):
  a = []
  while len(left) > 0 or len(right) > 0:
    if len(left) > 0 and len(right) > 0:
      if do_compare(left[0],right[0],args):
        a.append(left.pop(0))
      else:
        a.append(right.pop(0))
    elif len(left) > 0:
      a.append(left.pop(0))
    elif len(right) > 0:
      a.append(right.pop(0))
  return a


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

