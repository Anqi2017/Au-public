#!/usr/bin/python
import sys, os, inspect, argparse, multiprocessing, random
import PSLBasics, RangeBasics

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
  parser.add_argument('-o','--output',help="Write to output file, default is STDIN")
  parser.add_argument('--noheader',action='store_true')
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
  if args.threads > 1:
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
        if args.threads > 1:
          pool.apply_async(process_buffer,[buffer],callback=print_result)
        else:
          res = process_buffer(buffer)
          print_result(res)
      buffer = PSLBasics.MultiplePSLAlignments()
      if args.minimum_coverage > 1:
        buffer.set_minimum_coverage(args.minimum_coverage)
    last_name = e['qName']
    buffer.add_entry(PSLBasics.PSL(line.rstrip()))
  inf.close()
  if buffer.get_alignment_count() > 0:
    if args.threads > 1:
      pool.apply_async(process_buffer,[buffer],callback=print_result) # if we still have something left to do
    else:
      res = process_buffer(buffer)
      print_result(res)
  if args.threads > 1:
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

def print_result(results):
  if not results: return
  global lock
  global of
  lock.acquire()
  for result in results:
    of.write(result+"\n")
  lock.release()

def process_buffer(mpa):
  mpa.populate_query()
  #separate entires by Locus
  loci = RangeBasics.Loci()
  loci.set_use_direction(True)
  minimum_locus_distance = 400000
  for entry in mpa.entries:
    rng = RangeBasics.Bed(entry.value('tName'),entry.value('tStart')-minimum_locus_distance,entry.value('tEnd')+minimum_locus_distance,entry.value('strand'))
    rng.set_payload(entry)
    locus = RangeBasics.Locus()
    locus.set_use_direction(True)
    locus.add_member(rng)
    loci.add_locus(locus)
  loci.update_loci()
  outputs = []
  for locus in loci.loci:
    mpsl = PSLBasics.MultiplePSLAlignments()
    mpsl.set_minimum_coverage(20)
    for member in locus.members:
      mpsl.add_entry(member.get_payload())
    mpsl.populate_query()
    bac = mpsl.best_query()
    #if len(bac.segments) > 2:
    segtrimmed = bac.get_trimmed_entries()
    stitched = PSLBasics.stitch_query_trimmed_psl_entries(segtrimmed)
    #bac.print_report()
    #for segpsl in bac.segment_trimmed_entries:
    #  print str(segpsl.value('qStart'))+"\t"+str(segpsl.value('qEnd'))+"\t"+str(segpsl.value('tStart'))+"\t"+str(segpsl.value('tEnd'))
    ##print stitched.get_line()
    if not stitched.validate():
      sys.exit()
    outputs.append(stitched.get_line())
  return outputs
if __name__=="__main__":
  main()
