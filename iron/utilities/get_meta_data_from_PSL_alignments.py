#!/usr/bin/python
import sys, os, inspect, argparse, multiprocessing, random
import PSLBasics

##### Import local modules ####
#pythonfolder_loc = "IDP/pythonlib"
#cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
#if cmd_subfolder not in sys.path:
#  sys.path.insert(0,cmd_subfolder)
#import FusionPSLBasics
##### Done importing local modules ####

# The purpose of this script is to tear through PSL alignments and get
# details regarding whether they meet some criteria for being fusions or if they are multiply mapped etc..

# Global
lock = multiprocessing.Lock()
of = sys.stdout

def main():
  parser = argparse.ArgumentParser(description="Analyze ORDERED psl alignments of long reads.")
  parser.add_argument('psl_file',help="Alignment file. Must be ordered by query name. use - for stdin")
  parser.add_argument('--output',help="Write to output file, default is STDIN")
  parser.add_argument('--noheader',action='store_true')
  #parser.add_argument('--best',action='store_true')
  #parser.add_argument('--split',action='store_true')
  parser.add_argument('--minimum_coverage',type=int,help="Only consider alignments with at least this many bp aligned")
  parser.add_argument('--threads',type=int,default=multiprocessing.cpu_count(),help="INT default cpu_count")
  parser.add_argument('--tempbuffer',help="DIRECTORY store the results in a temporary file until they are ready to output.  suggest using /tmp if you don't know what to use")
  args = parser.parse_args()
  seen_names = set()
  last_name = ''
  buffer = PSLBasics.MultiplePSLAlignments()
  inf = sys.stdin
  if args.psl_file != '-':
    inf = open(args.psl_file)
  global of
  tname = None
  if args.tempbuffer:
    if not args.output:
      sys.stderr.write("ERROR if you want to buffer outputs in a temp file you need to specify a final output file.\n")
      sys.exit()
    rnum = random.randint(1,1000000000);
    tname = args.tempbuffer.rstrip('/')+'/weirathe.'+str(rnum)+'.meta'
    of = open(tname,'w')
  if args.output and not args.tempbuffer:
    of = open(args.output,'w')
  global lock
  if not args.noheader:
    lock.acquire()
    of.write("QueryName\tSegmentCount\tLocusCount\tHasOverlapped\tHasMultiplyMapped\n")
    lock.release()
  pool = multiprocessing.Pool(args.threads)
  for line in inf:
    e = PSLBasics.line_to_entry(line.rstrip())
    if e['qName'] != last_name: # we have a new name
      if e['qName'] in seen_names:
        sys.stderr.write("ERROR psl entries are not ordered by query name.\n")
        sys.exit()
      seen_names.add(e['qName'])
      if buffer.get_alignment_count() > 0:
        #process_buffer(buffer)
        pool.apply_async(process_buffer,[buffer],callback=print_result)
      buffer = PSLBasics.MultiplePSLAlignments()
      if args.minimum_coverage > 1:
        buffer.set_minimum_coverage(args.minimum_coverage)
    last_name = e['qName']
    buffer.add_entry(e)
  inf.close()
  if buffer.get_alignment_count() > 0:
    #process_buffer(buffer) # if we still have something left to do
    pool.apply_async(process_buffer,[buffer],callback=print_result) # if we still have something left to do
  pool.close()
  pool.join()
  of.close()
  if args.tempbuffer:
    of = open(args.output,'w')
    with open(tname) as inf:
      for line in inf:
        of.write(line)
    of.close()
    os.remove(tname)

def print_result(result):
  global lock
  global of
  lock.acquire()
  of.write(result+"\n")
  lock.release()

def process_buffer(mpa):
  mpa.populate_query()
  bac = mpa.best_query()
  has_overlapped = 0
  if bac.has_overlapped_segments(): has_overlapped = 1
  has_multiply_mapped = 0
  if bac.has_multiply_mapped_segments(): has_multiply_mapped = 1
  return bac.qName +"\t" + str(len(bac.segments)) + "\t" + str(bac.locus_count())+ "\t" + str(has_overlapped) + "\t" + str(has_multiply_mapped)

if __name__=="__main__":
  main()
