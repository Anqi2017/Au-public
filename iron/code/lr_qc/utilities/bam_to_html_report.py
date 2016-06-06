#!/usr/bin/python
import argparse, sys, os, time, re
from shutil import rmtree, copy, copytree
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

def main():
  #do our inputs
  args = do_inputs()

  if not args.output and not args.portable_output:
    sys.stderr.write("ERROR: must specify some kind of output\n")
    sys.exit()
  #Make sure rscript is installed
  try:
    cmd = 'Rscript --version'
    prscript = Popen(cmd.split(),stdout=PIPE,stderr=PIPE)
    rline = prscript.communicate()
    sys.stderr.write("Using Rscript version:\n")
    sys.stderr.write(rline[1].rstrip()+"\n")
  except:
    sys.stderr.write("ERROR: Rscript not installed\n")
    sys.exit()

  #Make sure python is installed
  try:
    cmd = 'python --version'
    prscript = Popen(cmd.split(),stdout=PIPE,stderr=PIPE)
    rline = prscript.communicate()
    sys.stderr.write("Using Python version:\n")
    sys.stderr.write(rline[1].rstrip()+"\n")
  except:
    sys.stderr.write("ERROR: python not installed\n")
    sys.exit()

  ## Check and see if directory for outputs exists
  if args.output:
    if os.path.isdir(args.output):
      sys.stderr.write("ERROR: output directory already exists.  Remove it to write to this location.\n")
      sys.exit()

  

  if not os.path.exists(args.tempdir+'/plots'):
    os.makedirs(args.tempdir+'/plots')
  if not os.path.exists(args.tempdir+'/data'):
    os.makedirs(args.tempdir+'/data')
  if not os.path.exists(args.tempdir+'/logs'):
    os.makedirs(args.tempdir+'/logs')
  # Make the plots for the page
  make_plots(args)
  make_html(args)
  of = open(args.tempdir+'/data/params.txt','w')
  for arg in vars(args):
    of.write(arg+"\t"+str(getattr(args,arg))+"\n")
  of.close()
  udir = os.path.dirname(os.path.realpath(__file__))
  if args.output:
    copytree(args.tempdir,args.output)
    cmd = 'python '+udir+'/make_solo_html.py '+args.output+'/report.html'
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    with open(args.output+'/portable_report.html','w') as of:
      for line in p.stdout:
        of.write(line)
    p.communicate()
  if args.portable_output:
    cmd = 'python '+udir+'/make_solo_html.py '+args.tempdir+'/report.html'
    sys.stderr.write(cmd+"\n")
    p = Popen(cmd.split(),stdout=PIPE)
    with open(args.portable_output,'w') as of:
      for line in p.stdout:
        of.write(line)
    p.communicate()
  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def make_plots(args):
  # Make plots
  udir = os.path.dirname(os.path.realpath(__file__))
  p = Pool(processes=args.threads)
  cmd = 'python '+udir+'/bam_to_gapped_alignment_plot.py '+args.input+' -r '+args.reference+' -o '+args.tempdir+'/plots/alignments_plot.png '+args.tempdir+'/plots/alignments_plot.pdf --output_stats '+args.tempdir+'/data/alignment_stats.txt'
  if args.max_query_overlap:
    cmd += ' --max_query_overlap '+str(args.max_query_overlap)
  if args.max_target_overlap:
    cmd += ' --max_target_overlap '+str(args.max_target_overlap)
  if args.max_query_gap:
    cmd += ' --max_query_gap '+str(args.max_query_gap)
  if args.max_target_gap:
    cmd += ' --max_target_gap '+str(args.max_target_gap)
  if args.required_fractional_improvement:
    cmd += ' --required_fractional_improvement '+str(args.required_fractional_improvement)
  sys.stderr.write("Making gapped alignment plot\n")
  sys.stderr.write(cmd+"\n")
  p.apply_async(mycall,args=(cmd,args.tempdir+'/logs/alignment',))  
  cmd = 'python '+udir+'/bam_to_context_error_plot.py '+args.input+' -r '+args.reference+' --target --output_raw '+args.tempdir+'/data/context_error_data.txt -o '+args.tempdir+'/plots/context_plot.png '+args.tempdir+'/plots/context_plot.pdf'
  if args.context_error_scale:
    cmd += ' --scale '+' '.join([str(x) for x in args.context_error_scale])
  if args.context_error_stopping_point:
    cmd += ' --stopping_point '+str(args.context_error_stopping_point)
  sys.stderr.write("Making context plot\n")
  sys.stderr.write(cmd+"\n")
  p.apply_async(mycall,args=(cmd,args.tempdir+'/logs/context_error',))  
  #call(cmd.split())  
  cmd = 'python '+udir+'/bam_to_alignment_error_plot.py '+args.input+' -r '+args.reference+' --output_stats '+args.tempdir+'/data/error_stats.txt --output_raw '+args.tempdir+'/data/error_data.txt -o '+args.tempdir+'/plots/alignment_error_plot.png '+args.tempdir+'/plots/alignment_error_plot.pdf'
  if args.alignment_error_scale:
    cmd += ' --scale '+' '.join([str(x) for x in args.alignment_error_scale])
  if args.alignment_error_max_length:
    cmd += ' --max_length '+str(args.alignment_error_max_length)
  sys.stderr.write("Making alignment error plot\n")
  sys.stderr.write(cmd+"\n")
  p.apply_async(mycall,args=(cmd,args.tempdir+'/logs/alignment_error',))  
  #call(cmd.split())  
  p.close()
  p.join()

def mycall(cmd,lfile):
  ofe = open(lfile+'.err','w')
  ofo = open(lfile+'.out','w')
  p = Popen(cmd.split(),stderr=ofe,stdout=ofo)
  p.communicate()
  ofe.close()
  ofo.close()
  return

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Create an output report",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
  parser.add_argument('-o','--output',help="OUTPUT Folder or STDOUT if not set")
  parser.add_argument('--portable_output',help="OUTPUT file in a portable html format")
  parser.add_argument('-r','--reference',help="Reference Fasta",required=True)
  parser.add_argument('--threads',type=int,default=1,help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")

  ### Parameters for alignment plots
  parser.add_argument('--max_query_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  parser.add_argument('--max_target_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  parser.add_argument('--max_query_gap',type=int,help="for testing gapped alignment advantge")
  parser.add_argument('--max_target_gap',type=int,default=500000,help="for testing gapped alignment advantage")
  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="require gapped alignment to be this much better (in alignment length) than single alignment to consider it.")
  
  ### Params for alignment error plot
  parser.add_argument('--alignment_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  parser.add_argument('--alignment_error_max_length',type=int,default=1000000,help="The maximum number of alignment bases to calculate error from")
  
  ### Params for context error plot
  parser.add_argument('--context_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  parser.add_argument('--context_error_stopping_point',type=int,default=10000,help="Sample at least this number of each context")
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

def make_html(args):
  #read in our alignment data
  mydate = time.strftime("%Y-%m-%d")
  a = {}
  with open(args.tempdir+'/data/alignment_stats.txt') as inf:
    for line in inf:
      (name,numstr)=line.rstrip().split("\t")
      a[name]=int(numstr)
  #read in our error data
  e = {}
  with open(args.tempdir+'/data/error_stats.txt') as inf:
    for line in inf:
      (name,numstr)=line.rstrip().split("\t")
      e[name]=int(numstr)

  #make our css directory
  if not os.path.exists(args.tempdir+'/css'):
    os.makedirs(args.tempdir+'/css')
  udir = os.path.dirname(os.path.realpath(__file__))
  #copy css into that directory
  copy(udir+'/../data/mystyle.css',args.tempdir+'/css/mystyle.css')
  of = open(args.tempdir+'/report.html','w')
  ostr = '''
<!DOCTYPE html>
<html>
<head>
<link rel="stylesheet" type="text/css" href="css/mystyle.css">
<title>Long Read Alignment and Error Report</title>
</head>
<body>
<div id="topline">
  <div id="date_block">
    <div>
    Generated on:
    </div>
    <div class="input_value" id="date">'''
  of.write(ostr)
  of.write(mydate)
  ostr = '''</div>
  </div>
  <div id="param_block">
    <div>Execution parmeters:</div>
    <div class="input_value">
    <a href="data/params.txt">params.txt</a>
    </div>
  </div>
</div>
<div class="clear"></div>
<div>
Long read alignment and error report for:
</div>
<div class="input_value" id="filename">'''
  of.write(ostr+"\n")
  of.write(args.input)
  ostr = '''</div>  
<hr>
<div id="alignment_analysis" class="subject_title"><table><tr><td class="c1">Alignment analysis</td><td class="c2"><span class="highlight">'''
  of.write(ostr)
  reads_aligned = perc(a['ALIGNED_READS'],a['TOTAL_READS'],1)
  of.write(reads_aligned)
  ostr = '''</span></td><td class="c3"><span class="highlight2">reads aligned</span></td><td class="c4"><span class="highlight">'''
  of.write(ostr)
  bases_aligned = perc(a['ALIGNED_BASES'],a['TOTAL_BASES'],1)
  of.write(bases_aligned)
  ostr = '''</span></td><td class="c5"><span class="highlight2">bases aligned <i>(of aligned reads)</i></span></td></tr></table>
</div>
<div id="alignment_block">
  <div id="alignment_stats">
    <table class="data_table">
      <tr class="rhead"><td colspan="3">Read Stats</td></tr>'''
  of.write(ostr+"\n")
  total_read_string = '<tr><td>Total reads</td><td>'+str(a['TOTAL_READS'])+'</td></td><td></td></tr>'
  of.write(total_read_string+"\n")
  unaligned_read_string = '<tr><td>- Unaligned reads</td><td>'+str(a['UNALIGNED_READS'])+'</td></td><td>'+perc(a['UNALIGNED_READS'],a['TOTAL_READS'],1)+'</td></tr>'
  of.write(unaligned_read_string+"\n")
  aligned_read_string = '<tr><td>- Aligned reads</td><td>'+str(a['ALIGNED_READS'])+'</td></td><td>'+perc(a['ALIGNED_READS'],a['TOTAL_READS'],1)+'</td></tr>'
  of.write(aligned_read_string+"\n")
  single_align_read_string = '<tr><td>--- Single-align reads</td><td>'+str(a['SINGLE_ALIGN_READS'])+'</td></td><td>'+perc(a['SINGLE_ALIGN_READS'],a['TOTAL_READS'],1)+'</td></tr>'
  of.write(single_align_read_string+"\n")
  gapped_align_read_string = '<tr><td>--- Gapped-align reads</td><td>'+str(a['GAPPED_ALIGN_READS'])+'</td></td><td>'+perc(a['GAPPED_ALIGN_READS'],a['TOTAL_READS'],2)+'</td></tr>'
  of.write(gapped_align_read_string+"\n")
  ostr='''<tr class="rhead"><td colspan="3">Base Stats <i>(among aligned reads)</i></td></tr>'''
  of.write(ostr+"\n")
  total_bases_string = '<tr><td>Total bases</td><td>'+str(a['TOTAL_BASES'])+'</td></td><td></td></tr>'
  of.write(total_bases_string+"\n")
  unaligned_bases_string = '<tr><td>- Unaligned bases</td><td>'+str(a['UNALIGNED_BASES'])+'</td><td>'+perc(a['UNALIGNED_BASES'],a['TOTAL_BASES'],1)+'</td></tr>'
  of.write(unaligned_bases_string+"\n")
  aligned_bases_string = '<tr><td>- Aligned bases</td><td>'+str(a['ALIGNED_BASES'])+'</td><td>'+perc(a['ALIGNED_BASES'],a['TOTAL_BASES'],1)+'</td></tr>'
  of.write(aligned_bases_string+"\n")
  single_align_bases_string = '<tr><td>--- Single-aligned bases</td><td>'+str(a['SINGLE_ALIGN_BASES'])+'</td><td>'+perc(a['SINGLE_ALIGN_BASES'],a['TOTAL_BASES'],1)+'</td></tr>'
  of.write(single_align_bases_string+"\n")
  gapped_align_bases_string = '<tr><td>--- Additional-aligned bases</td><td>'+str(a['GAPPED_ALIGN_BASES'])+'</td><td>'+perc(a['GAPPED_ALIGN_BASES'],a['TOTAL_BASES'],2)+'</td></tr>'
  of.write(gapped_align_bases_string+"\n")
  #gapped_improvement_string = '<tr><td>----- Gapped-improvement</td><td>'+"80000000"+'</td><td>'+"80%"+'</td></tr>'
  #of.write(gapped_improvement_string+"\n")
  ostr = '''</table>
    <div id="legend_block">
      <table id="align_legend">
        <tr><td>Unaligned</td><td id="unaligned_leg"></td></tr>
        <tr><td>Gapped alignment</td><td id="gapped_leg"></td></tr>
        <tr><td>Single alignment</td><td id="single_leg"></td></tr>
      </table>
    </div>
  </div>
  <div id="gapped_image_block">
    <div class="rhead">Summary [<a href="plots/alignments_plot.pdf">pdf</a>]</div>
    <img id="gapped_image" src="plots/alignments_plot.png">
  </div>   
</div>

<div class="clear"></div>
<hr>
<div class="subject_title">Error pattern analysis &nbsp;&nbsp;&nbsp;&nbsp;<span class="highlight">'''
  of.write(ostr+"\n")
  error_rate = perc(e['ANY_ERROR'],e['ALIGNMENT_BASES'],3)
  of.write(error_rate)
  ostr='''</span> <span class="highlight2">error rate</span></div>
<div class="subject_subtitle">&nbsp; &nbsp; &nbsp; based on aligned segments</div>
<div id="error_block">
  <div id="context_error_block">
    <div class="rhead">Error rates, given a target sequence [<a href="plots/context_plot.pdf">pdf</a>]</div>
    <img id="context_image" class="square_image" src="plots/context_plot.png">
  </div>
  <div class="clear"></div>
  <div id="error_stats">
    <table class="data_table">
      <tr class="rhead"><td colspan="3">Alignment stats</td></tr>'''
  of.write(ostr+"\n")
  best_alignments_sampled_string = '<tr><td>Best alignments sampled</td><td>'+str(e['ALIGNMENT_COUNT'])+'</td><td></td></tr>'
  of.write(best_alignments_sampled_string+"\n")
  ostr = '''<tr class="rhead"><td colspan="3">Base stats</td></tr>'''
  of.write(ostr+"\n")
  bases_analyzed_string = '<tr><td>Bases analyzed</td><td>'+str(e['ALIGNMENT_BASES'])+'</td><td></td></tr>'
  of.write(bases_analyzed_string+"\n")
  correctly_aligned_string = '<tr><td>- Correctly aligned bases</td><td>'+str(e['ALIGNMENT_BASES']-e['ANY_ERROR'])+'</td><td>'+perc((e['ALIGNMENT_BASES']-e['ANY_ERROR']),e['ALIGNMENT_BASES'],1)+'</td></tr>'
  of.write(correctly_aligned_string+"\n")
  total_error_string = '<tr><td>- Total error bases</td><td>'+str(e['ANY_ERROR'])+'</td><td>'+perc(e['ANY_ERROR'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(total_error_string+"\n")
  mismatched_string = '<tr><td>--- Mismatched bases</td><td>'+str(e['MISMATCHES'])+'</td><td>'+perc(e['MISMATCHES'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(mismatched_string+"\n")
  deletion_string = '<tr><td>--- Deletion bases</td><td>'+str(e['ANY_DELETION'])+'</td><td>'+perc(e['ANY_DELETION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(deletion_string+"\n")
  complete_deletion_string = '<tr><td>----- Complete deletion bases</td><td>'+str(e['COMPLETE_DELETION'])+'</td><td>'+perc(e['COMPLETE_DELETION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(complete_deletion_string+"\n")
  homopolymer_deletion_string = '<tr><td>----- Homopolymer deletion bases</td><td>'+str(e['HOMOPOLYMER_DELETION'])+'</td><td>'+perc(e['HOMOPOLYMER_DELETION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(homopolymer_deletion_string+"\n")
  insertion_string = '<tr><td>--- Insertion bases</td><td>'+str(e['ANY_INSERTION'])+'</td><td>'+perc(e['ANY_INSERTION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(insertion_string+"\n")
  complete_insertion_string = '<tr><td>----- Complete insertion bases</td><td>'+str(e['COMPLETE_INSERTION'])+'</td><td>'+perc(e['COMPLETE_INSERTION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(complete_insertion_string+"\n")
  homopolymer_insertion_string = '<tr><td>----- Homopolymer insertion bases</td><td>'+str(e['HOMOPOLYMER_INSERTION'])+'</td><td>'+perc(e['HOMOPOLYMER_INSERTION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(homopolymer_insertion_string+"\n")
  ostr = '''</table>
  </div>
  <div id="alignment_error_block">
    <div class="rhead">Alignment-based error rates [<a href="plots/alignment_error_plot.pdf">pdf<a/>]</div>
    <img class="square_image" src="plots/alignment_error_plot.png">
  </div>
</div>
<div class="clear"></div>
<hr>
<table class="header_table">
  <tr><td class="rhead" colspan="4">Raw data</td></tr>
  <tr>
    <td>Alignments stats raw report:</td>
    <td>Alignment errors data:</td>
    <td>Alignment error report:</td>
    <td>Contextual errors data:</td>
  </tr>
  <tr class="raw_files">
    <td><a href="data/alignment_stats.txt">alignment_stats.txt</a></td>
    <td><a href="data/error_data.txt">error_data.txt</a></td>
    <td><a href="data/error_stats.txt">error_stats.txt</a></td>
    <td><a href="data/context_error_data.txt">context_error_data.txt</a></td>
  </tr>
</table>
</body>
</html>
  '''
  of.write(ostr)

#Pre: numerator and denominator
#Post: percentage string
def perc(num,den,decimals=0):
  s = "{0:."+str(decimals)+"f}%"
  return s.format(100*float(num)/float(den))

if __name__=="__main__":
  main()
