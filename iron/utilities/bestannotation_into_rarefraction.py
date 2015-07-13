#!/usr/bin/python
import argparse, sys, re, multiprocessing, json, zlib
from random import shuffle

def main():
  parser = argparse.ArgumentParser(description='Convert a bestoutput from annotated_psl_with_gpd.py into a rarefraction curve.')
  parser.add_argument('--original_read_count',type=int,help="INT total number of reads we originally had.  If not specified, the curve will be relative to reads we have here.\n")
  parser.add_argument('--length_filter',help="INT,INT i.e. 0,1000 only count annotations with transcript lengths in this range 0 and 1000 included.")
  parser.add_argument('--min_read_count',type=int,default=1,help='INT minimum number of reads to consider the annotation')
  parser.add_argument('--threads',type=int,help='INT default = CPU count')
  parser.add_argument('--filter_names',help='FILENAME  newline delimited list of annotations.  Limit rarefraction curve to this subset of possible annotations.')
  parser.add_argument('--gpd_name',help='The gpd_name to create the rarefraction curve of.  Not this is a field in that bestannotation, and not necessarily a file name.  If you dont know what to put youll be prompted with options.')
  parser.add_argument('--full_only',action='store_true',help='Only base the output on full alignments')
  parser.add_argument('--annotation_type',default='gene',choices=['gene','transcript'],help='Make a curve by "gene" or "transcript"')
  parser.add_argument('bestannotation_input',help='FILENAME created by annotate_psl_with_gpd with the "--bestoutput" setting specificed.')
  args = parser.parse_args()
  cpus = multiprocessing.cpu_count()
  size_filter = False
  if args.length_filter:
    size_filter = [int(x) for x in args.length_filter.split(",")]
  if args.threads:
    cpus = args.threads    
  avail_names = get_avail_gpd_names(args.bestannotation_input)
  if not args.gpd_name or args.gpd_name not in avail_names:
    sys.stderr.write("Set a gpd name '--gpd_name' from one of the following:\n")
    for e in sorted(avail_names):
      sys.stderr.write("  "+e+"\n")
    return
  # parse the input file
  read_annotations = {}
  with open(args.bestannotation_input) as inf:
    lnum = 0
    for line in inf:
      lnum += 1
      if lnum == 1 and is_header(line): continue # its a header
      f = line.rstrip().split("\t")
      if args.full_only and f[9] != 'Full': continue
      rname = f[1]
      annot = f[6]
      if args.annotation_type == 'transcript':
        annot = f[7]
      read_annotations[rname] = {}
      read_annotations[rname]['annotation'] = annot
      read_annotations[rname]['size'] = int(f[15])

  # Lets clear out things not meeting our size filter
  if size_filter:  # if we are filtering on size do it here
    read_names = read_annotations.keys()
    for read in read_names:
      if read_annotations[read]['size'] < size_filter[0] or read_annotations[read]['size'] > size_filter[1]: 
        del read_annotations[read]

  # this is where we need to filter our reads based on a list of accepted annotations
  if args.filter_names:
    allowed_annots = set()
    with open(args.filter_names) as inf:
      for line in inf:
        allowed_annots.add(line.rstrip())
    reads = read_annotations.keys()
    for read in reads:
      annot = read_annotations[read]['annotation']
      if annot not in allowed_annots:
        del read_annotations[read]

  # After all filters we put back the original read count if we are using one.
  original_read_count = len(read_annotations)
  if args.original_read_count: 
    if args.original_read_count < len(read_annotations):
      sys.stderr.write("ERROR your original reads should not be less than your observed reads\n")
      return
    original_read_count = args.original_read_count

  num_to_add = original_read_count - len(read_annotations)
  for i in range(0,num_to_add):
    aname = 'arti_'+str(i)
    read_annotations[aname] = False

  # Now we are free to make our curve
  iter = 20000
  avg_samples = 20
  pool = multiprocessing.Pool(processes=cpus)
  check_vals = [0]
  next = 1
  sys.stderr.write(str(len(read_annotations))+ " total reads\n")
  while True:
    next = append_check(check_vals,next,len(read_annotations))
    if next > len(read_annotations):
      break
  check_vals += [len(read_annotations)]
  reads_string = zlib.compress(json.dumps(read_annotations))
  results = []
  job = 0
  job_total = len(check_vals)
  for size in check_vals:
    job += 1
    avg = pool.apply_async(get_level_average,[reads_string,size,avg_samples,args.min_read_count,job,job_total])
    results.append([size,avg])
    #get_level_average(reads_string,size,avg_samples,args.min_read_count)
  pool.close()
  pool.join()
  sys.stderr.write("\n")
  for entry in results:
    print str(entry[0])+"\t"+str(entry[1].get())

def append_check(clist,inval,max_size):
  current = inval
  next = inval*10
  working_next = next
  if next > max_size: working_next = max_size
  clist += range(current,working_next,current)
  return next

def get_level_average(reads_string,size,avg_samples,min_count,job,job_total):
  tot = 0
  num = 0
  reads = json.loads(zlib.decompress(reads_string))
  for i in range(0,avg_samples):
    num += 1
    annot_count = get_sample_count(reads,size,min_count)
    tot += annot_count
  avg = tot/num
  #print str(size) + "\t" + str(avg)
  sys.stderr.write("\rfinished ("+str(job)+'/'+str(job_total)+")")
  return avg

def get_sample_count(reads,size,min_count):
  read_names = reads.keys()
  shuffle(read_names)
  seen = {}
  if size > len(reads): size = len(reads)
  # get the annotations and how many times we've observed it
  for read in read_names[0:size]:
    if not reads[read]: continue # These are false entries 
    if reads[read]['annotation'] not in seen:
      seen[reads[read]['annotation']] = 0
    seen[reads[read]['annotation']] += 1
  # check that we hav reached the min_count
  satisfied = set()
  for annot in seen:
    if seen[annot] >= min_count:
      satisfied.add(annot)
  return len(satisfied)

def is_header(line):
  if re.match('^psl_entry_id',line):
    return True
  return False

def get_avail_gpd_names(infile):
  avail = set()
  with open(infile) as inf:
    lnum = 0
    for line in inf:
      line = line.rstrip()
      lnum += 1
      if lnum == 1 and is_header(line):
        continue # its the header
      #print line
      f = line.split("\t")
      avail.add(f[5])
  return avail

main()
