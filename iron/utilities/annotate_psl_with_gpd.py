#!/usr/bin/python
import sys, re, random, os, multiprocessing, argparse
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
  parser.add_argument('pslfile',nargs=1,help='FILENAME PSL filename to annotate')
  parser.add_argument('gpdfile',nargs='+',help='FILENAME(S) genePred file(s) providing annotations')
  args = parser.parse_args()

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

  # write the gpd and bed from the psl file
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
  ostring = ostring.rstrip("\t")
  # make a job list
  sys.stderr.write("Entering multiprocessing annotations on "+str(cpus)+" cpus\n")
  for j in range(1,num_jobs+1):
    p.apply_async(execute_job,[tdir,j,geneprednames,overlap_fraction])

  p.close()
  p.join()
  of.write(ostring+"\n")
  for j in range(1,num_jobs+1):
    with open(tdir+"/annotated_full_match."+str(j)+".txt") as inf:
      for line in inf:
        of.write(line.rstrip()+"\n")
  of.close()
  rmtree(tdir)

def execute_job(tdir,j,geneprednames,overlap_fraction):
  of = open(tdir+'/unannotated_full_match.'+str(j)+'.txt','w')
  for i in range(1,len(geneprednames)+1):
    jobid = str(j)+"_"+str(i)
    pre_annotate(tdir,tdir+'/reference.'+str(i)+'.bed',tdir+'/partreads.'+str(j)+'.bed',of,jobid,overlap_fraction)
  of.close()
  annotate(tdir,j,len(geneprednames))
  return jobid

def pre_annotate(tdir,ref_file,obs_file,of,jobid,overlap_fraction):
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
      f = line.rstrip().split("\t")
      ref_exon_count = int(f[6])
      obs_exon_count = int(f[13])
      if ref_exon_count != obs_exon_count:
        continue
      ref_exon = f[0]+':'+f[1]+'-'+f[2]
      obs_exon = f[8]+':'+f[9]+'-'+f[10]
      ref_id = int(f[3])
      obs_id = int(f[11])
      obs_name = f[12]
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
      results[obs_id][ref_id]['matches'].add(str(ref_exon)+'_'+str(obs_exon))
      reflen = int(f[2])-int(f[1])
      obslen = int(f[10])-int(f[9])
      overlap = int(f[15])
      smallest = sorted([float(overlap)/float(reflen), float(overlap)/float(obslen)])[0]
      results[obs_id][ref_id]['overlaps'][int(f[1])] = smallest
      
  for obs_id in results:
    for ref_id in results[obs_id]:
      match_count = len(results[obs_id][ref_id]['matches'])
      if match_count != results[obs_id][ref_id]['ref_exon_count']:
        continue
      exstarts = sorted(results[obs_id][ref_id]['overlaps'].keys())
      if len(exstarts) > 2:
        innerexstarts = exstarts[2:len(exstarts)-1]
        minover = 1.0
        for exstart in innerexstarts:
          if results[obs_id][ref_id]['overlaps'][exstart] < minover:
            minover = results[obs_id][ref_id]['overlaps'][exstart]
        if minover < overlap_fraction[1]: continue
      if results[obs_id][ref_id]['ref_strand'] == '+':
        if results[obs_id][ref_id]['overlaps'][exstarts[0]] < overlap_fraction[0] and overlap_fraction[0] > 0:
          continue
        if results[obs_id][ref_id]['overlaps'][exstarts[len(exstarts)-1]] < overlap_fraction[2] and overlap_fraction[2] > 0:
          continue
      if results[obs_id][ref_id]['ref_strand'] == '-':
        if results[obs_id][ref_id]['overlaps'][exstarts[0]] < overlap_fraction[2] and overlap_fraction[2] > 0:
          continue
        if results[obs_id][ref_id]['overlaps'][exstarts[len(exstarts)-1]] < overlap_fraction[0] and overlap_fraction[0] > 0:
          continue
      of.write(str(obs_id) + "\t" + results[obs_id][ref_id]['read_name'] + "\t" \
               + str(results[obs_id][ref_id]['obs_exon_count']) + "\t" \
               + str(results[obs_id][ref_id]['ref_exon_count']) + "\t" \
               + str(ref_id) + "\n")


def annotate(tdir,partid,colcount):
  # make columns of the annotations
  d = {}
  with open(tdir+"/entries.txt") as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      ref_id = f[2]
      if ref_id not in d: d[ref_id] = {}
      d[ref_id]['column'] = int(f[0])
      d[ref_id]['gene'] = f[3]
      d[ref_id]['transcript'] = f[4]
  reads = {}
  with open(tdir+'/unannotated_full_match.'+str(partid)+'.txt') as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      if f[0] not in reads:
        reads[f[0]] = {}
        reads[f[0]]['read_name'] = f[1]
        reads[f[0]]['read_exon_count'] = int(f[2])
        reads[f[0]]['results'] = []
        for i in range(0,colcount):
          temp = {}
          temp['genes'] = set()
          temp['transcripts'] = set()
          reads[f[0]]['results'].append(temp)
      reads[f[0]]['results'][d[f[4]]['column']-1]['genes'].add(d[f[4]]['gene'])
      reads[f[0]]['results'][d[f[4]]['column']-1]['transcripts'].add(d[f[4]]['transcript'])
  of = open(tdir+"/annotated_full_match."+str(partid)+".txt",'w')
  for readid in [str(y) for y in sorted([int(x) for x in reads.keys()])]:
    ostring = readid + "\t" + reads[readid]['read_name'] + "\t" + str(reads[readid]['read_exon_count']) + "\t" 
    for i in range(0,colcount):
      ostring += ','.join(reads[readid]['results'][i]['genes']) + "\t"
      ostring += ','.join(reads[readid]['results'][i]['transcripts']) + "\t"
    ostring.rstrip("\t")
    of.write(ostring+"\n")
  of.close()

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
      line = line.rstrip()
      entry = GenePredBasics.line_to_entry(line)
      of_entries.write(str(column_number)+ "\t" + simplenames[column_number-1] + "\t" + str(entry_number) + "\t" + entry['gene_name'] + "\t" + entry['name']+"\n")
      for i in range(0,len(entry['exonStarts'])):
        of_ref.write(entry['chrom'] + "\t" + str(entry['exonStarts'][i]) + "\t" \
                   + str(entry['exonEnds'][i]) + "\t" + str(entry_number) + "\t" \
                   + entry['gene_name'] + "\t" \
                   + entry['name'] + "\t" + str(len(entry['exonStarts'])) + "\t" \
                   + entry['strand'] \
                   + "\n")
    gfr.close()
    of_ref.close()
  of_entries.close()


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
      line = line.rstrip()
      entry = GenePredBasics.line_to_entry(line)
      for i in range(0,len(entry['exonStarts'])):
        of_bed.write(entry['chrom'] + "\t" + str(entry['exonStarts'][i]) + "\t" \
                     + str(entry['exonEnds'][i]) + "\t" + entry['name']+"\t" \
                     + entry['gene_name'] + "\t" + str(len(entry['exonStarts'])) + "\t" + entry['strand'] + "\n")   
    oc.close()
    of_bed.close()
  return num_jobs

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
