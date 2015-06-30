#!/usr/bin/python
import argparse, sys, os, random, re
from shutil import rmtree, copytree

# This will operate a pipeline to go from 
# Pre:
#     pacbio_raw .. the path to the bax.h5 or bas.h5 you want to process
#                please keep in mind that in this requires the accordingly 
#                named xml file to be in the parent directory as the results
#                are stored by pacbio.
#     output directory must not already exist.
# Post:
#     output direcotry by default is in the current working directory
#     and called pre-IDP_output/ but can be anything by settting --output
#         

def main():
  parser = argparse.ArgumentParser(description="take step 1 and step 2 outputs and make inputs suitable for IDP",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--step_1_folder',required=True,help='FOLDERNAME output of step one')
  parser.add_argument('--step_2_folder',required=True,help='FOLDERNAME output of step two')
  parser.add_argument('--lsc_replacement_threshold',type=float,default=0.9,help='Replace corrected with full length corrected when they are this fraction of the full length')
  parser.add_argument('--tempdir',default='/tmp',help='FOLDERNAME location of temporary directory')
  parser.add_argument('--output',default='output_pre-IDP_step-3_final',help='FOLDERNAME of output must not already exist')
  parser.add_argument('--save_tempdir',help='DIRECTORYNAME name of a directory to be created that has the temporary folder in it.')
  args = parser.parse_args()

  if os.path.isdir(args.output):
    sys.stderr.write("ERROR "+args.output+" output folder must not already exist.\n")
    sys.exit()
  if args.save_tempdir:
    if os.path.isdir(args.save_tempdir):
      sys.stderr.write("ERROR "+args.save_tempdir+" folder must not already exist.\n")
      sys.exit()

  tdir = setup_temporary_directory(args)  
  sys.stderr.write("working in "+tdir+"\n")

  if not os.path.exists(tdir+'/output'):
    os.makedirs(tdir+'/output')

  of_log = open(tdir+'/output/LOG','w')
  of_log.write(str(sys.argv)+"\n")

  sys.stderr.write("Replace LSC corrected with full when lengths are similar.\n")
  [cortot,correp] = execute_replacement(tdir,args.step_2_folder,args.lsc_replacement_threshold)
  of_log.write("LSC replacement threshold: "+str(args.lsc_replacement_threshold)+"\n")
  sys.stderr.write("Replaced "+str(correp)+ " of "+ str(cortot)+" where corrected length was similar to corrected full length\n")
  sys.stderr.write("Compile the non-redundant set of corrected and uncorrected long reads\n")
  [zhq,zlq,zsub] = make_nonredundant(tdir+'/output/lr_nonredundant.fa',args.step_1_folder.rstrip('/')+'/ccs_hq/ccs_hq.fa',tdir+'/swapped_corrected.fa',args.step_1_folder.rstrip('/')+'/subreads/subreads.fa')
  sys.stderr.write(str(zhq)+" high quality ccs reads\n")
  sys.stderr.write(str(zlq)+" lsc corrected reads\n")
  sys.stderr.write(str(zsub)+" longest subreads\n")
  of_log.write(str(zhq)+" high quality ccs reads\n")
  of_log.write(str(zlq)+" lsc corrected reads\n")
  of_log.write(str(zsub)+" longest subreads\n")
  added = make_isoform(tdir+'/output/lr_for_isoforms.fa',tdir+'/output/lr_nonredundant.fa',tdir+'/not_swapped_full.fa')
  sys.stderr.write("Added "+str(cortot-correp)+" full length lsc corrected sequences for isoform prediction fasta\n")
  of_log.write("Added "+str(cortot-correp)+" full length lsc corrected sequences for isoform prediction fasta\n")
  of_log.close()
  copytree(tdir+'/output',args.output)
  if args.save_tempdir:
    copytree(tdir,args.save_tempdir)
  rmtree(tdir)  

def make_isoform(output_fasta,nr_fasta,not_swapped_fasta):
  of = open(output_fasta,'w')
  with open(nr_fasta) as inf:
    for line in inf: of.write(line)
  with open(not_swapped_fasta) as inf:
    for line in inf: of.write(line)
  of.close()
  return 

def make_nonredundant(output_fasta,hq_fasta,lq_corrected_fasta,subread_fasta):
  hq = read_fasta_into_hash(hq_fasta)
  lq = read_fasta_into_hash(lq_corrected_fasta)
  sub = read_fasta_into_hash(subread_fasta)
  seen_names = set()
  zhq = 0
  zlq = 0
  zsub = 0
  of = open(output_fasta,'w')
  for name in hq:
    short_name = get_short_name(name)
    seen_names.add(short_name)
    zhq += 1
    of.write(">"+name+"\n"+hq[name]+"\n")
  for name in lq:
    short_name = get_short_name(name)
    seen_names.add(short_name)
    zlq += 1
    of.write(">"+name+"\n"+lq[name]+"\n")
  keep = {}
  for name in sub:
    short_name = get_short_name(name)
    if short_name not in seen_names:
      if short_name not in keep:
        keep[short_name] = {}
        keep[short_name]['name'] = ''
        keep[short_name]['len'] = 0
        keep[short_name]['seq'] = ''
      if len(sub[name]) > keep[short_name]: # new longest
        keep[short_name]['name'] = name
        keep[short_name]['len'] = len(sub[name])
        keep[short_name]['seq'] = sub[name]
  for short_name in keep:
    zsub +=1
    of.write(">"+keep[short_name]['name']+"\n"+keep[short_name]['seq']+"\n")
  of.close()
  return [zhq,zlq,zsub]

def get_short_name(name):
  m = re.match('^([^\/]+\/\d+)',name)
  short_name = name
  if m:
    short_name = m.group(1)
  else:
    sys.stderr.write("ERROR strange gene name "+name+"\n")
    sys.exit()
  return short_name

def execute_replacement(tdir,lsc_dir,thresh):
  full = read_fasta_into_hash(lsc_dir.rstrip('/')+'/full_LR.fa')
  corrected = read_fasta_into_hash(lsc_dir.rstrip('/')+'/corrected_LR.fa')
  # put full back into corrected when lengths are similar
  of_not_swap_full = open(tdir+'/not_swapped_full.fa','w')
  of_swap_corrected = open(tdir+'/swapped_corrected.fa','w')
  z = 0
  zswap = 0
  for name in corrected:
    z += 1
    short_name = name
    m = re.match('^(.*)\|[\d\.]+$',name)
    if m: short_name = m.group(1) 
    if short_name not in full:
      sys.stderr.write("ERROR: " + name + " not in full")
      sys.exit()
    full_len = len(full[short_name])
    corrected_len = len(corrected[name])
    if len(corrected[name]) > len(full[short_name]):
      sys.stderr.write("WARNING: length of corrected greater than length of full")
    if len(full[short_name])*thresh <= len(corrected[name]):
      of_swap_corrected.write('>'+name+"\n"+full[short_name]+"\n")
      zswap+=1
    else:
      of_swap_corrected.write('>'+name+"\n"+corrected[name]+"\n")
      of_not_swap_full.write('>'+short_name+"\n"+full[short_name]+"\n")
  of_not_swap_full.close()
  of_swap_corrected.close()
  return [z, zswap]

def setup_temporary_directory(args):
  if not os.path.isdir(args.tempdir):
    sys.stderr.write("ERROR invalid temporary directory "+args.tempdir+"\n")
    sys.exit()
  tdir = args.tempdir.rstrip('/')+'/'+'preidp.'+str(random.randint(1,10000000))
  os.makedirs(tdir)
  if not os.path.isdir(tdir):
    sys.stderr.write("ERROR failed to make working temp directory "+args.tempdir+"\n")
    sys.exit()
  return tdir

# pre:       A fasta file name
# post:      A dictionary of sequences keyed by name
# modifies:  none

def read_fasta_into_hash(fasta_filename):
  seqs = {}
  with open(fasta_filename) as f:
    seq = ''
    head = ''
    prog = re.compile('^>(.+)')
    for line in f:
      line = line.rstrip()
      m = prog.match(line)
      if m:
        if seq != '':
          seqs[head] = seq
          seq = ''
        head = m.group(1)
      else:
         seq = seq + line
    if seq != '':
      seqs[head] = seq
  return seqs

if __name__=="__main__":
  main()
