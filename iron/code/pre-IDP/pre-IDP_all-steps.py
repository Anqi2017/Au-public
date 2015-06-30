#!/usr/bin/python
import argparse, sys, os, random, multiprocessing, subprocess, re
from shutil import rmtree, copyfile, copytree

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
  parser = argparse.ArgumentParser(description="take the pacbio raw data to the necessary input for IDP",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--smrtanalysis_path',required=True,help='PATH to smrtanalysis 2.3.0 source directory')
  parser.add_argument('--pacbio_raw',required=True,help='FILENAME .bax.h5 or .bas.h5 REQUIRED')
  parser.add_argument('--threads',type=int,default=0,help='INT number of threads to use')
  parser.add_argument('--tempdir',default='/tmp',help='FOLDERNAME location of temporary directory')
  parser.add_argument('--output',default='output_pre-IDP_all-steps',help='FOLDERNAME of output must not already exist')
  parser.add_argument('--ccs_hq_acc',type=int,default=95,help='INT accuracy of high quality ccs reads')
  parser.add_argument('--ccs_hq_passes',type=int,default=2,help='INT number of passes for high quality ccs reads')
  parser.add_argument('--ccs_lq_acc',type=int,default=75,help='INT accuracy of low quality ccs reads')
  parser.add_argument('--ccs_lq_passes',type=int,default=0,help='INT number of passes for low quality ccs reads')
  parser.add_argument('--subreads_acc',type=int,default=75,help='INT minimum accuracy of subreads')
  parser.add_argument('--save_tempdir',help='DIRECTORYNAME name of a directory to be created that has the temporary folder in it.')
  parser.add_argument('--short_reads',nargs='+',required=True,help='READFILE or READFILE_1 READFILE_2 for paired end')
  parser.add_argument('--short_read_type',required=True,choices=['FASTA','FASTQ'],help='FASTA or FASTQ')
  parser.add_argument('--lsc_replacement_threshold',type=float,default=0.9,help='Replace corrected with full length corrected when they are this fraction of the full length')
  args = parser.parse_args()

  if args.threads==0:
    args.threads = multiprocessing.cpu_count()

  if os.path.isdir(args.output):
    sys.stderr.write("ERROR "+args.output+" output folder must not already exist\n")
    sys.exit()

  tdir = setup_temporary_directory(args)  
  sys.stderr.write("working in "+tdir+"\n")

  if not os.path.exists(tdir+'/output'):
    os.makedirs(tdir+'/output')


  cmd1  = './pre-IDP_step-1_process_pacbio.py --smrtanalysis_path '+args.smrtanalysis_path+' '
  cmd1 += '--pacbio_raw '+args.pacbio_raw+' '
  cmd1 += '--output '+tdir+'/output/1_processed_pacbio '
  cmd1 += '--ccs_hq_acc '+str(args.ccs_hq_acc)+' '
  cmd1 += '--ccs_hq_passes '+str(args.ccs_hq_passes)+' '
  cmd1 += '--ccs_lq_acc '+str(args.ccs_lq_acc)+' '
  cmd1 += '--ccs_lq_passes '+str(args.ccs_lq_passes)+' '
  cmd1 += '--subreads_acc '+str(args.subreads_acc)+' '
  cmd1 += '--tempdir '+args.tempdir+' '
  cmd1 += '--threads '+str(args.threads)
  sys.stderr.write("STEP 1\n"+cmd1+"\n")
  subprocess.call(cmd1,shell=True)
  
  cmd2  = './pre-IDP_step-2_error_correction.py --step_1_folder '+tdir+'/output/1_processed_pacbio '
  cmd2 += '--short_reads '+' '.join([str(x) for x in args.short_reads])+' '
  cmd2 += '--short_read_type '+args.short_read_type+' '
  cmd2 += '--threads '+str(args.threads)+' '
  cmd2 += '--tempdir '+args.tempdir+' '
  cmd2 += '--output '+tdir+'/output/2_error_corrected '
  sys.stderr.write("STEP 2\n"+cmd2+"\n")
  subprocess.call(cmd2,shell=True)

  cmd3  = './pre-IDP_step-3_make_long_read_set.py --step_1_folder '+tdir+'/output/1_processed_pacbio '
  cmd3 += '--step_2_folder '+tdir+'/output/2_error_corrected '
  cmd3 += '--lsc_replacement_threshold '+str(args.lsc_replacement_threshold)+' '
  cmd3 += '--tempdir '+args.tempdir+' '
  cmd3 += '--output '+tdir+'/output/3_for_idp_long_reads'
  sys.stderr.write("STEP 3\n"+cmd3+"\n")
  subprocess.call(cmd3,shell=True)

  copytree(tdir+'/output',args.output)
  rmtree(tdir)

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

if __name__=="__main__":
  main()
