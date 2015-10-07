#!/usr/bin/python
import sys, re, random, os, multiprocessing
import argparse
import json
import GenePredBasics, PSLBasics, FileBasics, BigFileBasics
from shutil import rmtree

# Pre:  A long read psl file.  Any number of genepred files in the format "Gene Predictions and RefSeq Genes with Gene Names".
#       gzipped files are supported
# Post: A table with one row for each long read, and two columns for each annotation file
#       The columns for each annotation file correspond to Genes, and then Transcripts
#       It is possible for multiple genes or multiple transcripts to be reported back,
#       and in those cases they will be comma separated.
#       <read_id (unique)> <read_name> <read exon count> <gpd_1:gene> <gpd_1:transcript> ... <gpd_N:gene> < gpd_N:transcript>

def main():
  parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--output','-o',help='FILENAME output file, default: STDOUT')
  parser.add_argument('--rawoutput',help='FILENAME to write a db friendly output before any chosing best hits or reformating the report takes place')
  parser.add_argument('--bestoutput',help='FILENAME to write a db friendly output after chosing best hits but before reformating the report takes place')
  parser.add_argument('--gpdout',help='FILENAME location to output the psl files genePred conversion')
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default='/tmp',help='DIRECTORY temdir where a temporary directoy will be made and used')
  group.add_argument('--specific_tempdir',help='DIRECTORY temdir the exact tempdir to be used')
  parser.add_argument('--closegap',type=int,default=68,help='INT close gaps less than or equal to this')
  parser.add_argument('--jobsize',type=int,default=500000,help='not a very important parameter, it just says how many jobs to send to a thread at a time INT')
  parser.add_argument('--threads',type=int,help='INT number of threads default: cpu_count')
  parser.add_argument('--minoverlap',default='0,0.8,0',help='FLOATLIST First exon, inner exons, Last exon requirements.  Zero indicates that any overlap will suffice.')
  parser.add_argument('--mincoverage',type=float,default=0,help='FLOAT fraction of overall coverage we want to make a call')
  parser.add_argument('--minmatchbases',type=int,default=100,help='INT minimum number of bases needed to call a match default (100)')
  parser.add_argument('--debug',action='store_true',help='dont remove temporary files on execute when debugging')
  parser.add_argument('--input_is_gpd',action='store_true',help='instead of PSL the input is already genepred.')
  parser.add_argument('pslfile',help='FILENAME PSL filename to annotate - for STDIN')
  parser.add_argument('gpdfile',nargs='+',help='FILENAME(S) genePred file(s) providing annotations')
  args = parser.parse_args()

  # There is a lot going on so lets fill up this params dictionary with
  #   variables that we will be using.  Starting with the command line args
  params = {}
  params['args'] = args
  if args.specific_tempdir:
    params['tdir'] = args.specific_tempdir.rstrip('/')
  else:
    tbase = args.tempdir
    if args.tempdir:  
      tbase = args.tempdir.rstrip('/')
      if not os.path.exists(tbase):
        sys.stderr.write("Error: temp directory does not exist. "+args.tempdir)
        return
    # Store our actual temporary directory that we will create and write to
    params['tdir'] = tbase + '/weirathe.ea' + str(random.randint(1,1000000000))
  # Convert the overlap fraction string to an array of floats
  params['overlap_fraction'] = [float(x) for x in args.minoverlap.split(',')]
  sys.stderr.write("Temp directory: "+params['tdir']+"\n")
  # make the temporary directory
  if not os.path.exists(params['tdir']): 
    os.makedirs(params['tdir'])

  # where to output results to
  of = sys.stdout
  if args.output:
    of = open(args.output,'w')

  # where to output results to
  if args.gpdout:
    #just making sure we can write here so we don't get a surprise later
    gpdof = open(args.gpdout,'w')

  sys.stderr.write("pslfile: "+args.pslfile+"\n")
  sys.stderr.write("Converting psl file to gpd\n")

  # write the gpd from the psl file
  if not args.input_is_gpd:
    print params['args'].closegap
    parse_pslfile(params['tdir'],params['args'].pslfile,params['args'].closegap)
  else:
    parse_gpdfile(params['tdir'],params['args'].pslfile,params['args'].closegap)
  # save the genepred if we want it
  if args.gpdout:
    sys.stderr.write("saving genepred conversion of psl file to: "+args.gpdout+"\n")
    with open(params['tdir']+'/longreads.gpd') as gf:
      for line in gf: gpdof.write(line)
    gpdof.close()  

  # break the new gpd into jobs
  sys.stderr.write("Splitting job\n")
  params['num_jobs'] = break_gpdfile(params['tdir'],params['args'].jobsize)

  simplenames = [] #names to be used to label columns
  for file in params['args'].gpdfile:
    m = re.search('([^\/]+$)',file)
    name = file
    sys.stderr.write("  "+name+"\n")
    if m: name = m.group(1)
    simplenames.append(name)
  params['simplenames'] = simplenames

  # convert reference genepreds to bed fies
  sys.stderr.write("Parsing reference file\n")
  parse_refgpd(params['tdir'],params['args'].gpdfile,params['simplenames'])

  if not params['args'].threads: params['args'].threads = multiprocessing.cpu_count()
  if params['args'].threads > 1:
    p = multiprocessing.Pool(processes=params['args'].threads)
  # make a job list
  sys.stderr.write("Entering multiprocessing annotations "+str(params['num_jobs'])+" jobs on "+str(params['args'].threads)+" cpus\n")
  for j in range(1,params['num_jobs']+1):
    # The business happens here with execute job.
    if params['args'].threads > 1:
      p.apply_async(execute_job,[params['tdir'],j,params['args'].gpdfile,params['overlap_fraction'],params['args'].minmatchbases])
    else:
      execute_job(params['tdir'],j,params['args'].gpdfile,params['overlap_fraction'],params['args'].minmatchbases)
    #execute_job(params['tdir'],j,params['args'].gpdfile,params['overlap_fraction'],params['args'].minmatchbases)
  if params['args'].threads > 1:
    p.close()
    p.join()

  
  # Print out the raw data here if we want it
  #   and save our best match per read/gpd 
  read_gpd = {}
  columns = {}
  ostring = "psl_entry_id\tread_name\tread_exons\treference_exons\tgpd_column_number\t"
  ostring += "gpd_name\tgene_name\ttranscript_name\tfragment_report\talignment_classification\tsplit_alignment\ttotal_aligned_bases\ttotal_aligned_exons\tlongest_fragment_bases\tlongest_fragment_exons\treference_length\n"
  if params['args'].rawoutput:
    of_raw = open(params['args'].rawoutput,'w')
    of_raw.write(ostring)
  for j in range(1,params['num_jobs']+1):
    with open(params['tdir']+"/annotated_match."+str(j)+".txt") as inf:
      for line in inf:
        f = line.rstrip("\n").split("\t")
        my_gpd = f[5]
        my_read = f[1]
        columns[int(f[4])] = my_gpd
        if my_read not in read_gpd:  
          read_gpd[my_read] = {}
        if my_gpd not in read_gpd[my_read]:
          read_gpd[my_read][my_gpd] = {}
          read_gpd[my_read][my_gpd]['Full'] = {}
          read_gpd[my_read][my_gpd]['Full']['best_hit'] = False
          read_gpd[my_read][my_gpd]['Full']['matches'] = 0
          read_gpd[my_read][my_gpd]['Partial'] = {}
          read_gpd[my_read][my_gpd]['Partial']['best_hit'] = False
          read_gpd[my_read][my_gpd]['Partial']['matches'] = 0
          read_gpd[my_read][my_gpd]['Best'] = False
        total_matches = int(f[11])
        if f[9] == 'Full' and total_matches > read_gpd[my_read][my_gpd]['Full']['matches']:
          read_gpd[my_read][my_gpd]['Full']['matches'] = total_matches
          read_gpd[my_read][my_gpd]['Full']['best_hit'] = f
        if f[9] == 'Partial' and total_matches > read_gpd[my_read][my_gpd]['Partial']['matches']:
          read_gpd[my_read][my_gpd]['Partial']['matches'] = total_matches
          read_gpd[my_read][my_gpd]['Partial']['best_hit'] = f
        if params['args'].rawoutput: of_raw.write(line.rstrip("\n")+"\n")

  if params['args'].bestoutput: 
    ostring = "psl_entry_id\tread_name\tread_exons\treference_exons\tgpd_column_number\t"
    ostring += "gpd_name\tgene_name\ttranscript_name\tfragment_report\talignment_classification\tsplit_alignment\ttotal_aligned_bases\ttotal_aligned_exons\tlongest_fragment_bases\tlongest_fragment_exons\treference_length\n"
    ofbest = open(params['args'].bestoutput,'w')
    ofbest.write(ostring)
  for read in read_gpd:
    for gpd in read_gpd[read]:
      if read_gpd[read][gpd]['Full']['matches'] > 0:
        cov = get_cov(read_gpd[read][gpd]['Full']['best_hit'][11],read_gpd[read][gpd]['Full']['best_hit'][15])
        if not params['args'].mincoverage or cov >= params['args'].mincoverage:
          if params['args'].bestoutput: ofbest.write("\t".join(read_gpd[read][gpd]['Full']['best_hit'])+"\n") 
          read_gpd[read][gpd]['Best'] = read_gpd[read][gpd]['Full']['best_hit']
      elif read_gpd[read][gpd]['Partial']['matches'] > 0:
        cov = get_cov(read_gpd[read][gpd]['Partial']['best_hit'][11],read_gpd[read][gpd]['Partial']['best_hit'][15])
        if not params['args'].mincoverage or cov >= params['args'].mincoverage:
          if params['args'].bestoutput: ofbest.write("\t".join(read_gpd[read][gpd]['Partial']['best_hit'])+"\n")
          read_gpd[read][gpd]['Best'] = read_gpd[read][gpd]['Partial']['best_hit']
  if params['args'].bestoutput: ofbest.close()
  #Now lets do the final report form output
  
  colnums = sorted(columns.keys())
  ostring = "read\t"
  for colnum in colnums:
    ostring += 'genes:'+columns[colnum]+"\t"
    ostring += 'transcripts:'+columns[colnum]+"\t"
    ostring += 'classification:'+columns[colnum]+"\t"
  ostring = ostring[:-1]
  of.write(ostring+"\n")
  for read in read_gpd:
    ostring = read + "\t"
    seen = 0
    for colnum in colnums:
      done = 0
      if columns[colnum] in read_gpd[read]:
        if read_gpd[read][columns[colnum]]['Best']:
          ostring += read_gpd[read][columns[colnum]]['Best'][6] + "\t"
          ostring += read_gpd[read][columns[colnum]]['Best'][7] + "\t"
          ostring += read_gpd[read][columns[colnum]]['Best'][9] + "\t"
          done = 1
          seen = 1
      if done == 0:
        ostring += "\t"
    ostring = ostring[:-1]
    if seen == 1:
      of.write(ostring+"\n")
  of.close()
  if not params['args'].debug and not params['args'].specific_tempdir: rmtree(params['tdir'])


def get_cov(f1,f2):
  if f1 <= 0 or f2 <= 0: return 0
  return min(float(f1)/float(f2),float(f2)/float(f1))

# This is how we call the process of working on one of our results
# Pre: Temporary Directory, job number, list of genepred files, overlap_fraction
#      where overlap fraction is an array of the required overlap for the
#      first, internal, and last exons
# 
def execute_job(tdir,j,geneprednames,overlap_fraction,min_match_bp):
  of = open(tdir+'/annotated_match.'+str(j)+'.txt','w')
  for i in range(1,len(geneprednames)+1):
    # Assign jobid as the job number, and the genepred column number
    jobid = str(j)+"_"+str(i)
    pre_annotate(tdir,tdir+'/reference.'+str(i)+'.bed',tdir+'/partreads.'+str(j)+'.bed',of,jobid,overlap_fraction,min_match_bp)
  of.close()
  #annotate(tdir,j,geneprednames)
  return jobid

# This pre-annotate is where we actually overlap the files
#   We do an intersection opertation, then we read through the
#   intersection seeing if it meets criteria for matching
# Pre:  Temporary Directory
#       Reference genepred bed file
#       long reads job bed file
#       output file handle for 'unannotated full match'
#       jobid jobnumber underscore column number (genepred file number)
#       overlap fraction
# Post: writes intersection of beds for each bed file
#       into unannotated_match.(job).txt
#   1.  Reads PSL entry number
#   2.  Read name
#   3.  Observed exon count
#   4.  Reference exon count
#   5.  Reference entry number
#   6.  Consecutive exon match(s)

def pre_annotate(tdir,ref_file,obs_file,of,jobid,overlap_fraction,min_match_bp):
  #print ref_file + "\t" + obs_file
  cmd = "bedtools intersect -wo " 
  #if overlap_fraction > 0: cmd += "-r -f "+str(overlap_fraction)
  cmd += " -a " + ref_file + " -b " + obs_file + " > " + tdir + "/intersect."+jobid+".bed" 
  #print cmd
  os.system(cmd)
  # now parse the intersection file

  # check for a full length match
  results = {}
  with open(tdir+'/intersect.'+jobid+'.bed') as inf:
    for line in inf:
      f = line.rstrip("\n").split("\t")
      ref_exon_count = int(f[6])
      obs_exon_count = int(f[14])
      #if ref_exon_count != obs_exon_count:
      #  continue
      ref_exon = f[0]+':'+f[1]+'-'+f[2]
      obs_exon = f[9]+':'+f[10]+'-'+f[11]
      ref_id = int(f[3])
      obs_id = int(f[12])
      obs_name = f[13]
      if obs_id not in results:
        results[obs_id] = {}
      if ref_id not in results[obs_id]:
        results[obs_id][ref_id] = {}
        results[obs_id][ref_id]['read_name'] = obs_name
        results[obs_id][ref_id]['matches'] = set()
        results[obs_id][ref_id]['ref_exon_count'] = ref_exon_count
        results[obs_id][ref_id]['obs_exon_count'] = obs_exon_count
        results[obs_id][ref_id]['name'] = f[4]
        results[obs_id][ref_id]['transcript'] = f[5]
        results[obs_id][ref_id]['ref_strand'] = f[7]
        results[obs_id][ref_id]['exons_ref'] = {}
        results[obs_id][ref_id]['exon_overlap'] = {}
      results[obs_id][ref_id]['matches'].add(str(ref_exon)+'_'+str(obs_exon))
      reflen = int(f[2])-int(f[1])
      obslen = int(f[11])-int(f[10])
      overlap = int(f[17])
      # get the overlap fraction
      if reflen == 0 or obslen == 0:
        sys.stderr.write(line+"\n")
        smallest = 0
      else:
        smallest = sorted([float(overlap)/float(reflen), float(overlap)/float(obslen)])[0]
      ref_exon_number = int(f[8])
      obs_exon_number = int(f[16])
      results[obs_id][ref_id]['exons_ref'][ref_exon_number] = obs_exon_number
      if ref_exon_number not in results[obs_id][ref_id]['exon_overlap']:
        results[obs_id][ref_id]['exon_overlap'][ref_exon_number] = {}
      results[obs_id][ref_id]['exon_overlap'][ref_exon_number][obs_exon_number] = {}
      results[obs_id][ref_id]['exon_overlap'][ref_exon_number][obs_exon_number]['bp'] = overlap
      results[obs_id][ref_id]['exon_overlap'][ref_exon_number][obs_exon_number]['frac'] = smallest


  d = {}
  with open(tdir+"/entries.txt") as inf:
    for line in inf:
      f = line.rstrip("\n").split("\t")
      ref_id = int(f[2])
      if ref_id not in d: d[ref_id] = {}
      d[ref_id]['column'] = int(f[0])
      d[ref_id]['gpdname'] = f[1]
      d[ref_id]['gene'] = f[3]
      d[ref_id]['transcript'] = f[4]
      d[ref_id]['length'] = f[5]

  #Go through the results and find the best consecutive exons
  for obs_id in results:
    for ref_id in results[obs_id]:
      overlap_data = results[obs_id][ref_id]['exon_overlap']
      refnums = sorted(results[obs_id][ref_id]['exons_ref'].keys())
      prevrefval = refnums[0]
      prevobsval = results[obs_id][ref_id]['exons_ref'][refnums[0]]
      best = []
      allconsec = []
      best.append([prevrefval,prevobsval])
      for n in refnums[1:]:
        obs = results[obs_id][ref_id]['exons_ref'][n]
        if n != prevrefval+1 or obs != prevobsval+1:
          allconsec.append(best)
          best = []
        best.append([n,obs])
        prevrefval = n
        prevobsval = obs
      if len(best) > 0: 
        allconsec.append(best)
      # now the allconsec contains all the consecutive bests
      passing_consec = {}
      match_bases = 0
      for consec in allconsec:
        fracs = [overlap_data[x[0]][x[1]]['frac'] for x in consec]
        totalbps = 0
        for bp in [overlap_data[x[0]][x[1]]['bp'] for x in consec]: totalbps += bp
        passing = True
        if len(consec) > 2:
          for frac in fracs[1:len(fracs)-1]:
            if frac < overlap_fraction[1]:
              passing = False
        if results[obs_id][ref_id]['ref_strand'] == '+':
          if fracs[0] < overlap_fraction[0] and overlap_fraction[0] > 0:
            passing = False
          if fracs[len(fracs)-1] < overlap_fraction[2] and overlap_fraction[2] > 0:
            passing = False
        if results[obs_id][ref_id]['ref_strand'] == '-':
          if fracs[0] < overlap_fraction[2] and overlap_fraction[2] > 0:
            passing = False
          if fracs[len(fracs)-1] < overlap_fraction[0] and overlap_fraction[0] > 0:
            passing = False
        if passing:
          passing_consec[json.dumps(consec)] = {}
          passing_consec[json.dumps(consec)]['bp'] = totalbps
          match_bases += totalbps
          passing_consec[json.dumps(consec)]['exons'] = len(consec)

      total_aligned_bases = 0
      total_aligned_exons = 0
      longest_fragment_aligned_bases = 0
      longest_fragment_exon_count = 0
      for v in passing_consec:
        total_aligned_bases += passing_consec[v]['bp']
        total_aligned_exons += passing_consec[v]['exons']
        # consider longest fragment by base pairs for now
        #   the alternative would be exon count
        if passing_consec[v]['bp'] > longest_fragment_aligned_bases:
          longest_fragment_aligned_bases = passing_consec[v]['bp']
          longest_fragment_exon_count = passing_consec[v]['exons']
      if len(passing_consec) == 0: continue #make sure we passed our criteria      
      if match_bases < min_match_bp: continue
      matchstring =  ",".join([str(passing_consec[x]['exons'])+":"+str(passing_consec[x]['bp']) for x in passing_consec])
      matchtype = 'Full'
      if results[obs_id][ref_id]['ref_exon_count'] != results[obs_id][ref_id]['obs_exon_count'] and results[obs_id][ref_id]['ref_exon_count'] != len(passing_consec):
        matchtype = 'Partial'
      gappedtype = 'N'
      if len(passing_consec) > 1:
        gappedtype = 'Y'
      of.write(str(obs_id) + "\t" + results[obs_id][ref_id]['read_name'] + "\t" \
               + str(results[obs_id][ref_id]['obs_exon_count']) + "\t" \
               + str(results[obs_id][ref_id]['ref_exon_count']) + "\t" \
               + str(d[ref_id]['column']) + "\t" + d[ref_id]['gpdname'] + "\t" \
               + d[ref_id]['gene'] + "\t" + d[ref_id]['transcript'] + "\t" \
               + matchstring + "\t"  \
               + matchtype + "\t" + gappedtype + "\t" \
               + str(total_aligned_bases) + "\t" + str(total_aligned_exons) + "\t" \
               + str(longest_fragment_aligned_bases) + "\t" + str(longest_fragment_exon_count) + "\t" \
               + str(d[ref_id]['length']) + "\n")

#Parse the reference genepred(s) into bed file
#Pre: temporary directory, genepredfilenames, simplenames
#     temporary directory - path to temporary directory
#     genepredfilenames - list of reference gpd filenames
#     simplenames - list of genepred short names
#Post: Bed file with the following format
#  1.  chrom
#  2.  start 0-base
#  3.  end 1-base
#  4.  reference gpd entry line number
#  5.  gene name
#  6.  transcript name
#  7.  number of exons in reference gpd entry
#  8.  strand
#  9.  exon number
#  Writes to two places.  entries.txt and reference.(column_number).bed
def parse_refgpd(tdir,geneprednames,simplenames):
  # get the reference genepreds ready to use in work
  column_number = 0
  entry_number = 0
  of_entries = open(tdir+"/entries.txt",'w')
  for file in geneprednames:
    column_number += 1
    of_ref = open(tdir+"/reference."+str(column_number)+".bed",'w')
    gfr = FileBasics.GenericFileReader(file)
    while True:
      line = gfr.readline()
      if not line: break
      if re.match('^#',line): continue
      entry_number += 1
      line = line.rstrip("\n")
      entry = GenePredBasics.line_to_entry(line)
      entry_length = 0
      for i in range(0,len(entry['exonStarts'])): entry_length += entry['exonEnds'][i]-entry['exonStarts'][i]
      of_entries.write(str(column_number)+ "\t" + simplenames[column_number-1] + "\t" + str(entry_number) + "\t" + entry['gene_name'] + "\t" + entry['name']+"\t"+str(entry_length)+"\n")
      exon_number = 0
      for i in range(0,len(entry['exonStarts'])):
        exon_number += 1
        of_ref.write(entry['chrom'] + "\t" + str(entry['exonStarts'][i]) + "\t" \
                   + str(entry['exonEnds'][i]) + "\t" + str(entry_number) + "\t" \
                   + entry['gene_name'] + "\t" \
                   + entry['name'] + "\t" + str(len(entry['exonStarts'])) + "\t" \
                   + entry['strand'] + "\t" + str(exon_number) \
                   + "\n")
    gfr.close()
    of_ref.close()
  of_entries.close()

# Break the genpred into jobs
# Pre:  Temporary directory, job size (int)
# Post:  Write a bed file from each jobs segement of the genepred file
#        Bed file format is as follows:
#    1.  chrom
#    2.  start
#    3.  end
#    4.  PSL entry number
#    5.  read name
#    6.  number of exons
#    7.  strand
#    8.  exon_number
#  Write the bed files into partreads.(job).bed
#  job is an integer 1-based
#  we return the number of jobs also
def break_gpdfile(tdir,job_size):
  bfcr = BigFileBasics.BigFileChunkReader(tdir+'/longreads.gpd')
  bfcr.set_chunk_size_bytes(job_size)
  num_jobs = bfcr.chunk_count
  for i in range(0,bfcr.chunk_count):
    oc = bfcr.open_chunk(i)
    job = i+1
    of_bed = open(tdir+'/partreads.'+str(job)+'.bed','w')
    while True:
      line = oc.read_line()
      if not line: break
      line = line.rstrip("\n")
      entry = GenePredBasics.line_to_entry(line)
      exon_number = 0
      for i in range(0,len(entry['exonStarts'])):
        exon_number += 1
        of_bed.write(entry['chrom'] + "\t" + str(entry['exonStarts'][i]) + "\t" \
                     + str(entry['exonEnds'][i]) + "\t" + entry['name']+"\t" \
                     + entry['gene_name'] + "\t" + str(len(entry['exonStarts'])) + "\t" \
                     + entry['strand'] + "\t" + str(exon_number) + "\n")   
    oc.close()
    of_bed.close()
  return num_jobs

#Write the genepred
# Pre: temp directory, the psl file, smoothing factor (min intron size)
# Post: into longreads.gpd we write the genepred line
def parse_pslfile(tdir,pslfile,smoothing_factor):
  # Go through the long reads and make a genepred
  if pslfile != '-':
    fr = FileBasics.GenericFileReader(pslfile)
  else:
    fr = sys.stdin
  seennames = {}
  longreadnumber = 0
  of_gpd = open(tdir+'/longreads.gpd','w')
  while True:
    line = fr.readline()
    if not line: break
    if re.match('^#',line): #skip comments
      continue
    longreadnumber += 1
    gpd_line = PSLBasics.convert_entry_to_genepred_line(PSLBasics.line_to_entry(line.rstrip()))
    if not gpd_line:
      sys.stderr.write("Warning: malformed psl for "+readname+"\n")
      continue
    entry = GenePredBasics.smooth_gaps( \
              GenePredBasics.line_to_entry(gpd_line),smoothing_factor)
    readname = entry['name']
    if readname in seennames:
      sys.stderr.write("Warning: repeat name '"+readname+"'\n")
    #set our first name to our bin
    entry['name'] = str(longreadnumber)
    gline = GenePredBasics.entry_to_line(entry)
    of_gpd.write(gline+"\n")
  fr.close()
  of_gpd.close()

#Write the genepred
# Pre: temp directory, the psl file, smoothing factor (min intron size)
# Post: into longreads.gpd we write the genepred line
def parse_gpdfile(tdir,gpdfile,smoothing_factor):
  # Go through the long reads and make a genepred
  if gpdfile != '-':
    fr = FileBasics.GenericFileReader(gpdfile)
  else:
    fr = sys.stdin
  seennames = {}
  longreadnumber = 0
  of_gpd = open(tdir+'/longreads.gpd','w')
  while True:
    line = fr.readline()
    if not line: break
    if re.match('^#',line): #skip comments
      continue
    longreadnumber += 1
    entry = GenePredBasics.smooth_gaps( \
              GenePredBasics.line_to_entry(line.rstrip()) \
              ,smoothing_factor)
    readname = entry['name']
    if readname in seennames:
      sys.stderr.write("Warning: repeat name '"+readname+"'\n")
    #set our first name to our bin
    entry['name'] = str(longreadnumber)
    gline = GenePredBasics.entry_to_line(entry)
    of_gpd.write(gline+"\n")
  fr.close()
  of_gpd.close()



main()
