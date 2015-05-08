#!/usr/bin/python
import sys, re, random, os, multiprocessing, argparse, json
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
  parser = argparse.ArgumentParser()
  parser.add_argument('-o',nargs='?',help='FILENAME output file, default: STDOUT')
  parser.add_argument('--gpdout',nargs='?',help='FILENAME location to output the psl files genePred conversion')
  parser.add_argument('--tempdir',nargs='?',help='DIRECTORY temdir, default: /tmp')
  parser.add_argument('--closegap',nargs='?',help='INT close gaps less than or equal to this, default: 68')
  parser.add_argument('--jobsize',nargs='?',help='INT default: 5000000')
  parser.add_argument('--minoverlap',nargs='?',help='FLOATLIST First exon, inner exons, Last exon requirements.  Zero indicates that any overlap will suffice. default: 0,0.8,0')
  parser.add_argument('--minmatchbases',nargs='?',help='INT minimum number of bases needed to call a match default (100)')
  parser.add_argument('--debug',action='store_true')
  parser.add_argument('pslfile',nargs=1,help='FILENAME PSL filename to annotate')
  parser.add_argument('gpdfile',nargs='+',help='FILENAME(S) genePred file(s) providing annotations')
  args = parser.parse_args()

  min_match_bp = 100
  if args.minmatchbases:
    min_match_bp = int(args.minmatchbp)
  #temp_dir_base = '/localscratch/weirathe'
  temp_dir_base = '/tmp'
  if args.tempdir:  
    temp_dir_base = args.tempdir.rstrip('/')
    if not os.path.exists(temp_dir_base):
      sys.stderr.write("Error: temp directory does not exist. "+temp_dir_base)
      return
  smoothing_factor = 68  # distance to close gaps
  if args.closegap:  smoothing_factor = int(args.closegap)

  job_size = 5000000
  if args.jobsize: job_size = int(args.jobsize)

  inner_overlap_fraction = 0.8  # reciprocial overlap required except for 
                                # the last exon of multiexon genes
  first_overlap_fraction = 0
  last_overlap_fraction = 0  # single exon genes are subject to the smaller of the first or lap settings
  if args.minoverlap: [first_overlap_fraction, inner_overlap_fraction, last_overlap_fraction] = [float(x) for x in args.minoverlap.split(',')]
  overlap_fraction = [first_overlap_fraction, inner_overlap_fraction, last_overlap_fraction]

  temp_name = 'weirathe.ea' + str(random.randint(1,1000000000))
  tdir = temp_dir_base + '/' + temp_name
  sys.stderr.write("Temp directory: "+tdir+"\n")

  # where to output results to
  of = sys.stdout
  if args.o:
    of = open(args.o,'w')

  # make the temporary directory
  if not os.path.exists(tdir): 
    os.makedirs(tdir)

  # where to output results to
  if args.gpdout:
    #just making sure we can write here so we don't get a surprise later
    gpdof = open(args.gpdout,'w')

  pslfile = args.pslfile[0]
  sys.stderr.write("pslfile: "+pslfile+"\n")
  geneprednames = args.gpdfile
  sys.stderr.write("Converting psl file to gpd\n")

  # write the gpd from the psl file
  parse_pslfile(tdir,pslfile,smoothing_factor)

  # save the genepred if we want it
  if args.gpdout:
    sys.stderr.write("saving genepred conversion of psl file to: "+args.gpdout+"\n")
    with open(tdir+'/longreads.gpd') as gf:
      for line in gf: gpdof.write(line)
    gpdof.close()  

  # break the new gpd into jobs
  sys.stderr.write("Splitting job\n")
  num_jobs = break_gpdfile(tdir,job_size)

  simplenames = [] #names to be used to label columns
  for file in geneprednames:
    m = re.search('([^\/]+$)',file)
    name = file
    sys.stderr.write("  "+name+"\n")
    if m: name = m.group(1)
    simplenames.append(name)

  # convert reference genepreds to bed fies
  sys.stderr.write("Parsing reference file\n")
  parse_refgpd(tdir,geneprednames,simplenames)

  cpus = multiprocessing.cpu_count()
  p = multiprocessing.Pool(processes=cpus)
  ostring = "read_id\tread_name\tread_exons\t"
  for name in simplenames:
    ostring += name+":genes\t"+name+":transcripts\t"
  ostring = ostring[:-1]
  # make a job list
  sys.stderr.write("Entering multiprocessing annotations "+str(num_jobs)+" jobs on "+str(cpus)+" cpus\n")
  for j in range(1,num_jobs+1):
    # The business happens here with execute job.
    #p.apply_async(execute_job,[tdir,j,geneprednames,overlap_fraction,min_match_bp])
    execute_job(tdir,j,geneprednames,overlap_fraction,min_match_bp)
  p.close()
  p.join()
  of.write(ostring+"\n")
  for j in range(1,num_jobs+1):
    with open(tdir+"/annotated_full_match."+str(j)+".txt") as inf:
      for line in inf:
        of.write(line.rstrip("\n")+"\n")
  of.close()
  if not args.debug: rmtree(tdir)

# This is how we call the process of working on one of our results
# Pre: Temporary Directory, job number, list of genepred files, overlap_fraction
#      where overlap fraction is an array of the required overlap for the
#      first, internal, and last exons
# 
def execute_job(tdir,j,geneprednames,overlap_fraction,min_match_bp):
  of = open(tdir+'/unannotated_match.'+str(j)+'.txt','w')
  for i in range(1,len(geneprednames)+1):
    # Assign jobid as the job number, and the genepred column number
    jobid = str(j)+"_"+str(i)
    pre_annotate(tdir,tdir+'/reference.'+str(i)+'.bed',tdir+'/partreads.'+str(j)+'.bed',of,jobid,overlap_fraction,min_match_bp)
  of.close()
  annotate(tdir,j,len(geneprednames))
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
        results[obs_id][ref_id]['overlaps'] = {}
        results[obs_id][ref_id]['ref_strand'] = f[7]
        results[obs_id][ref_id]['exons_ref'] = {}
        results[obs_id][ref_id]['exon_overlap'] = {}
        results[obs_id][ref_id]['exon_overlap'] = {}
      results[obs_id][ref_id]['matches'].add(str(ref_exon)+'_'+str(obs_exon))
      reflen = int(f[2])-int(f[1])
      obslen = int(f[11])-int(f[10])
      overlap = int(f[17])
      smallest = sorted([float(overlap)/float(reflen), float(overlap)/float(obslen)])[0]
      results[obs_id][ref_id]['overlaps'][int(f[1])] = smallest
      ref_exon_number = int(f[8])
      obs_exon_number = int(f[16])
      results[obs_id][ref_id]['exons_ref'][ref_exon_number] = obs_exon_number
      if ref_exon_number not in results[obs_id][ref_id]['exon_overlap']:
        results[obs_id][ref_id]['exon_overlap'][ref_exon_number] = {}
      results[obs_id][ref_id]['exon_overlap'][ref_exon_number][obs_exon_number] = {}
      results[obs_id][ref_id]['exon_overlap'][ref_exon_number][obs_exon_number]['bp'] = overlap
      results[obs_id][ref_id]['exon_overlap'][ref_exon_number][obs_exon_number]['frac'] = smallest


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
               + str(ref_id) + "\t"+ matchstring + "\t"  \
               + matchtype + "\t" + gappedtype + "\n")

def annotate(tdir,partid,colcount):
  #print colcount
  # make columns of the annotations
  d = {}
  with open(tdir+"/entries.txt") as inf:
    for line in inf:
      f = line.rstrip("\n").split("\t")
      ref_id = f[2]
      if ref_id not in d: d[ref_id] = {}
      d[ref_id]['column'] = int(f[0])
      d[ref_id]['gene'] = f[3]
      d[ref_id]['transcript'] = f[4]
  reads = {}
  with open(tdir+'/unannotated_match.'+str(partid)+'.txt') as inf:
    for line in inf:
      f = line.rstrip("\n").split("\t")
      rid = f[0]
      if rid not in reads:
        reads[rid] = {}
        reads[rid]['read_name'] = f[1]
        reads[rid]['read_exon_count'] = int(f[2])
        #reads[f[0]]['ref_exon_count'] = int(f[3])
        #reads[f[0]]['matchstring'] = f[5]
        #reads[f[0]]['partial'] = f[6]
        #reads[f[0]]['gapped'] = f[7]
        #reads[f[0]]['ref_exon_count'] = int(f[3])
        reads[rid]['results'] = []
      # now reads all have a dataset
      columns = []
      for i in range(0,colcount):
        temp = {}
        #temp['genes'] = set()
        #temp['transcripts'] = set()
        columns.append(temp)
      columns[d[f[4]]['column']-1]['genes'] = d[f[4]]['gene']
      columns[d[f[4]]['column']-1]['transcripts'] = d[f[4]]['transcript']
      columns[d[f[4]]['column']-1]['ref_exon_count'] = int(f[3])
      columns[d[f[4]]['column']-1]['matchstring'] = f[5]
      columns[d[f[4]]['column']-1]['partial'] = f[6]
      columns[d[f[4]]['column']-1]['gapped'] = f[7]
      reads[rid]['results'].append(columns)
  of = open(tdir+"/annotated_full_match."+str(partid)+".txt",'w')
  for readid in [str(y) for y in sorted([int(x) for x in reads.keys()])]:
    for i in range(0,len(reads[readid]['results'])):
      ostring = readid + "\t" + reads[readid]['read_name'] + "\t" + str(reads[readid]['read_exon_count']) + "\t" 
      for j in range(0,colcount):
        ostring += reads[readid]['results'][i][j]['genes'] + "\t"
        ostring += reads[readid]['results'][i][j]['transcripts'] + "\t"
        ostring += str(reads[readid]['results'][i][j]['ref_exon_count']) + "\t"
        ostring += reads[readid]['results'][i][j]['matchstring'] + "\t"
        ostring += reads[readid]['results'][i][j]['partial'] + "\t"
        ostring += reads[readid]['results'][i][j]['gapped'] + "\t"
      ostring=ostring[:-1]
      of.write(ostring+"\n")
  of.close()

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
      of_entries.write(str(column_number)+ "\t" + simplenames[column_number-1] + "\t" + str(entry_number) + "\t" + entry['gene_name'] + "\t" + entry['name']+"\n")
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
  fr = FileBasics.GenericFileReader(pslfile)
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
              GenePredBasics.line_to_entry( \
                PSLBasics.convert_entry_to_genepred_line( \
                  PSLBasics.line_to_entry(line.rstrip() \
            ))),smoothing_factor)
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
