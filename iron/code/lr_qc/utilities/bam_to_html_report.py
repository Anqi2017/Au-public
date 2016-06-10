#!/usr/bin/python
import argparse, sys, os, time, re, gzip, locale
from shutil import rmtree, copy, copytree
from multiprocessing import cpu_count, Pool
from tempfile import mkdtemp, gettempdir
from subprocess import Popen, PIPE

# read count
rcnt = 0

locale.setlocale(locale.LC_ALL,'en_US')

def main():
  #do our inputs
  args = do_inputs()

  if not args.output and not args.portable_output:
    sys.stderr.write("ERROR: must specify some kind of output\n")
    sys.exit()

  if args.no_reference:
    sys.stderr.write("WARNING: No reference specified.  Will be unable to calcualte error pattern\n")

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

  ## Extract data that can be realized from the bam
  #make_data_bam(args)

  ## Extract data that can be realized from the bam and reference
  #if args.reference:
  #  make_data_bam_reference(args)

  ## Extract data that can be realized from bam and reference annotation
  #if args.annotation:
  #  make_data_bam_annotation(args)

  # Create the output HTML
  make_html(args)

  of = open(args.tempdir+'/data/params.txt','w')
  for arg in vars(args):
    of.write(arg+"\t"+str(getattr(args,arg))+"\n")
  of.close()
  udir = os.path.dirname(os.path.realpath(__file__))
  sys.exit()
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

def make_data_bam(args):
  # Get the data necessary for making tables and reports
  udir = os.path.dirname(os.path.realpath(__file__))
  cmd = 'python '+udir+'/bam_traversal.py '+args.input+' -o '+args.tempdir+'/data/ '
  cmd += ' --threads '+str(args.threads)+' '
  if args.min_aligned_bases:
    cmd += ' --min_aligned_bases '+str(args.min_aligned_bases)
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
  sys.stderr.write("Traverse bam for alignment analysis\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/bam_traversal')

  cmd = "gpd_to_bed_depth.py "+args.tempdir+'/data/best.sorted.gpd.gz -o '+args.tempdir+'/data/depth.sorted.bed.gz'
  sys.stderr.write("Generate the depth bed for the mapped reads\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/gpd_to_depth')

  cmd = 'python '+udir+"/gpd_loci_analysis.py "+args.tempdir+'/data/best.sorted.gpd.gz -o '+args.tempdir+'/data/loci-all.bed.gz --output_loci '+args.tempdir+'/data/loci.bed.gz'
  cmd += ' --threads '+str(args.threads)+' '
  if args.min_depth:
    cmd += ' --min_depth '+str(args.min_depth)
  if args.min_depth:
    cmd += ' --min_coverage_at_depth '+str(args.min_coverage_at_depth)
  if args.min_exon_count:
    cmd += ' --min_exon_count '+str(args.min_exon_count)
  sys.stderr.write("Approximate loci and mapped read distributions among them.\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/gpd_to_loci')
  global rcnt #read count
  rcnt = 0
  tinf = gzip.open(args.tempdir+'/data/lengths.txt.gz')
  for line in tinf:  rcnt += 1
  tinf.close()
  cmd = 'python '+udir+"/locus_bed_to_rarefraction.py "+args.tempdir+'/data/loci.bed.gz -o '+args.tempdir+'/data/locus_rarefraction.txt'
  cmd += ' --threads '+str(args.threads)+' '
  cmd += ' --original_read_count '+str(rcnt)+' '
  sys.stderr.write("Make rarefraction curve\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/loci_rarefraction')

  sys.stderr.write("Make locus rarefraction plot\n")
  for ext in ['png','pdf']:
    cmd = 'Rscript '+udir+'/plot_annotation_rarefractions.r '+\
             args.tempdir+'/plots/locus_rarefraction.'+ext+' '+\
             'locus'+' '+\
             args.tempdir+'/data/locus_rarefraction.txt '+\
             '#FF000088 '
    sys.stderr.write(cmd+"\n")
    mycall(cmd,args.tempdir+'/logs/plot_locus_rarefraction_'+ext)


  cmd = "python "+udir+'/make_alignment_plot.py '+args.tempdir+'/data/lengths.txt.gz '
  cmd += ' --output_stats '+args.tempdir+'/data/alignment_stats.txt '
  cmd += ' --output '+args.tempdir+'/plots/alignments.png '
  cmd += args.tempdir+'/plots/alignments.pdf'
  sys.stderr.write("Make alignment plots\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/alignment_plot')

  # Make depth reports
  sys.stderr.write("Making depth reports\n")
  cmd = "python "+udir+'/depth_to_coverage_report.py '+args.tempdir+'/data/depth.sorted.bed.gz '+args.tempdir+'/data/chrlens.txt -o '+args.tempdir+'/data'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/depth_to_coverage')

  # do the depth graphs
  sys.stderr.write("Making coverage plots\n")
  cmd = 'Rscript '+udir+'/plot_chr_depth.r  '+args.tempdir+'/data/line_plot_table.txt.gz '+args.tempdir+'/data/total_distro_table.txt.gz '+args.tempdir+'/data/chr_distro_table.txt.gz '+args.tempdir+'/plots/covgraph.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/covgraph_png')
  cmd = 'Rscript '+udir+'/plot_chr_depth.r  '+args.tempdir+'/data/line_plot_table.txt.gz '+args.tempdir+'/data/total_distro_table.txt.gz '+args.tempdir+'/data/chr_distro_table.txt.gz '+args.tempdir+'/plots/covgraph.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/covgraph_pdf')

  # do depth plots
  sys.stderr.write("Making chr depth plots\n")
  cmd = 'Rscript '+udir+'/plot_depthmap.r '+args.tempdir+'/data/depth.sorted.bed.gz '+args.tempdir+'/data/chrlens.txt '+args.tempdir+'/plots/perchrdepth.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/perchr_depth_png')
  cmd = 'Rscript '+udir+'/plot_depthmap.r '+args.tempdir+'/data/depth.sorted.bed.gz '+args.tempdir+'/data/chrlens.txt '+args.tempdir+'/plots/perchrdepth.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/perchr_depth_pdf')
  return  

def make_data_bam_reference(args):
  # make the context error plots
  udir = os.path.dirname(os.path.realpath(__file__))
  cmd = 'python '+udir+'/bam_to_context_error_plot.py '+args.input+' -r '+args.reference+' --target --output_raw '+args.tempdir+'/data/context_error_data.txt -o '+args.tempdir+'/plots/context_plot.png '+args.tempdir+'/plots/context_plot.pdf'
  if args.context_error_scale:
    cmd += ' --scale '+' '.join([str(x) for x in args.context_error_scale])
  if args.context_error_stopping_point:
    cmd += ' --stopping_point '+str(args.context_error_stopping_point)
  sys.stderr.write("Making context plot\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/context_error')

  cmd = 'python '+udir+'/bam_to_alignment_error_plot.py '+args.input+' -r '+args.reference+' --output_stats '+args.tempdir+'/data/error_stats.txt --output_raw '+args.tempdir+'/data/error_data.txt -o '+args.tempdir+'/plots/alignment_error_plot.png '+args.tempdir+'/plots/alignment_error_plot.pdf'
  if args.alignment_error_scale:
    cmd += ' --scale '+' '.join([str(x) for x in args.alignment_error_scale])
  if args.alignment_error_max_length:
    cmd += ' --max_length '+str(args.alignment_error_max_length)
  sys.stderr.write("Making alignment error plot\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/alignment_error')  

def make_data_bam_annotation(args):
  # make the context error plots
  udir = os.path.dirname(os.path.realpath(__file__))
  cmd = 'gpd_annotate.py '+args.tempdir+'/data/best.sorted.gpd.gz -r '+args.annotation+' -o '+args.tempdir+'/data/annotbest.txt.gz'
  if args.threads:
    cmd += ' --threads '+str(args.threads)
  sys.stderr.write("Annotating reads\n")
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/gpd_annotate')

  sys.stderr.write("Writing rarefraction curves\n")
  global rcnt
  cmd =  'python '+udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --gene -o '+args.tempdir+'/data/gene_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/gene_rarefraction')
  cmd =  'python '+udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --transcript -o '+args.tempdir+'/data/transcript_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/transcript_rarefraction')
  cmd =  'python '+udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --full --gene -o '+args.tempdir+'/data/gene_full_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/gene_full_rarefraction')
  cmd =  'python '+udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --full --transcript -o '+args.tempdir+'/data/transcript_full_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/transcript_full_rarefraction')

  # now make the plots
  for type in ['gene','transcript']:
    for ext in ['png','pdf']:
      cmd = 'Rscript '+udir+'/plot_annotation_rarefractions.r '+\
             args.tempdir+'/plots/'+type+'_rarefraction.'+ext+' '+\
             type+' '+\
             args.tempdir+'/data/'+type+'_rarefraction.txt '+\
             '#FF000088 '+\
              args.tempdir+'/data/'+type+'_full_rarefraction.txt '+\
             '#0000FF88 '
      sys.stderr.write(cmd+"\n")
      mycall(cmd,args.tempdir+'/logs/plot_'+type+'_rarefraction_'+ext)
  return


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
  group1 = parser.add_mutually_exclusive_group(required=True)
  group1.add_argument('-r','--reference',help="Reference Fasta")
  group1.add_argument('--no_reference',action='store_true',help="No Reference Fasta")
  parser.add_argument('--annotation',help="Reference annotation genePred")
  parser.add_argument('--threads',type=int,default=1,help="INT number of threads to run. Default is system cpu count")
  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")

  ### Parameters for alignment plots
  parser.add_argument('--min_aligned_bases',type=int,default=50,help="for analysizing alignment, minimum bases to consider")
  parser.add_argument('--max_query_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  parser.add_argument('--max_target_overlap',type=int,default=10,help="for testing gapped alignment advantage")
  parser.add_argument('--max_query_gap',type=int,help="for testing gapped alignment advantge")
  parser.add_argument('--max_target_gap',type=int,default=500000,help="for testing gapped alignment advantage")
  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="require gapped alignment to be this much better (in alignment length) than single alignment to consider it.")
  
  ### Parameters for locus analysis
  parser.add_argument('--min_depth',type=float,default=1.5,help="require this or more read depth to consider locus")
  parser.add_argument('--min_coverage_at_depth',type=float,default=0.8,help="require at leas this much of the read be covered at min_depth")
  parser.add_argument('--min_exon_count',type=int,default=2,help="Require at least this many exons in a read to consider assignment to a locus")

  ### Params for alignment error plot
  parser.add_argument('--alignment_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  parser.add_argument('--alignment_error_max_length',type=int,default=100000,help="The maximum number of alignment bases to calculate error from")
  
  ### Params for context error plot
  parser.add_argument('--context_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
  parser.add_argument('--context_error_stopping_point',type=int,default=1000,help="Sample at least this number of each context")
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

  # read in our coverage data
  totref = 0
  with open(args.tempdir+'/data/chrlens.txt') as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      totref+=int(f[1])
  inf = gzip.open(args.tempdir+'/data/depth.sorted.bed.gz')
  totcov = 0
  for line in inf:
    f = line.rstrip().split("\t")
    totcov += int(f[2])-int(f[1])
  inf.close()

  #get our locus count
  inf = gzip.open(args.tempdir+'/data/loci.bed.gz')
  locuscount = 0
  for line in inf:
    locuscount += 1
  inf.close()

  #get our annotation counts
  genefull = 0
  geneany = 0
  txfull = 0
  txany = 0
  inf = gzip.open(args.tempdir+'/data/annotbest.txt.gz')
  genes_f = {}
  genes_a = {}
  txs_f = {}
  txs_a = {}
  for line in inf:
    f = line.rstrip().split("\t")
    g = f[2]
    t = f[3]
    if g not in genes_a: genes_a[g] = 0
    genes_a[g]+=1
    if t not in txs_a: txs_a[t] = 0
    txs_a[t]+=1
    if f[4] == 'full':
      if g not in genes_f: genes_f[g] = 0
      genes_f[g]+=1
      if t not in txs_f: txs_f[t] = 0
      txs_f[t]+=1
  inf.close()
  genefull = len(genes_f.keys())
  geneany = len(genes_a.keys())
  txfull = len(txs_f.keys())
  txany = len(txs_a.keys())
  

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
  ostr = '''
    </div>
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
  ostr = '''
</div>  
<hr>
<div id="alignment_analysis" class="subject_title"><table><tr><td class="c1">Alignment analysis</td><td class="c2"><span class="highlight">'''
  of.write(ostr)
  reads_aligned = perc(a['ALIGNED_READS'],a['TOTAL_READS'],1)
  of.write(reads_aligned)
  ostr = '''
  </span></td><td class="c3"><span class="highlight2">reads aligned</span></td><td class="c4"><span class="highlight">'''
  of.write(ostr)
  bases_aligned = perc(a['ALIGNED_BASES'],a['TOTAL_BASES'],1)
  of.write(bases_aligned)
  ostr = '''
  </span></td><td class="c5"><span class="highlight2">bases aligned <i>(of aligned reads)</i></span></td></tr></table>
</div>
<div id="alignment_block">
  <div id="alignment_stats">
    <table class="data_table">
      <tr class="rhead"><td colspan="3">Read Stats</td></tr>'''
  of.write(ostr+"\n")
  total_read_string = '<tr><td>Total reads</td><td>'+str(addcommas(a['TOTAL_READS']))+'</td></td><td></td></tr>'
  of.write(total_read_string+"\n")
  unaligned_read_string = '<tr><td>- Unaligned reads</td><td>'+str(addcommas(a['UNALIGNED_READS']))+'</td></td><td>'+perc(a['UNALIGNED_READS'],a['TOTAL_READS'],1)+'</td></tr>'
  of.write(unaligned_read_string+"\n")
  aligned_read_string = '<tr><td>- Aligned reads</td><td>'+str(addcommas(a['ALIGNED_READS']))+'</td></td><td>'+perc(a['ALIGNED_READS'],a['TOTAL_READS'],1)+'</td></tr>'
  of.write(aligned_read_string+"\n")
  single_align_read_string = '<tr><td>--- Single-align reads</td><td>'+str(addcommas(a['SINGLE_ALIGN_READS']))+'</td></td><td>'+perc(a['SINGLE_ALIGN_READS'],a['TOTAL_READS'],1)+'</td></tr>'
  of.write(single_align_read_string+"\n")
  gapped_align_read_string = '<tr><td>--- Gapped-align reads</td><td>'+str(addcommas(a['GAPPED_ALIGN_READS']))+'</td></td><td>'+perc(a['GAPPED_ALIGN_READS'],a['TOTAL_READS'],2)+'</td></tr>'
  of.write(gapped_align_read_string+"\n")
  gapped_align_read_string = '<tr><td>--- Chimeric reads</td><td>'+str(addcommas(a['CHIMERA_ALIGN_READS']))+'</td></td><td>'+perc(a['CHIMERA_ALIGN_READS'],a['TOTAL_READS'],2)+'</td></tr>'
  of.write(gapped_align_read_string+"\n")
  gapped_align_read_string = '<tr><td>----- Trans-chimeric reads</td><td>'+str(addcommas(a['TRANSCHIMERA_ALIGN_READS']))+'</td></td><td>'+perc(a['TRANSCHIMERA_ALIGN_READS'],a['TOTAL_READS'],2)+'</td></tr>'
  of.write(gapped_align_read_string+"\n")
  gapped_align_read_string = '<tr><td>----- Self-chimeric reads</td><td>'+str(addcommas(a['SELFCHIMERA_ALIGN_READS']))+'</td></td><td>'+perc(a['SELFCHIMERA_ALIGN_READS'],a['TOTAL_READS'],2)+'</td></tr>'
  of.write(gapped_align_read_string+"\n")
  ostr='''
      <tr class="rhead"><td colspan="3">Base Stats <i>(of aligned reads)</i></td></tr>'''
  of.write(ostr+"\n")
  total_bases_string = '<tr><td>Total bases</td><td>'+str(addcommas(a['TOTAL_BASES']))+'</td></td><td></td></tr>'
  of.write(total_bases_string+"\n")
  unaligned_bases_string = '<tr><td>- Unaligned bases</td><td>'+str(addcommas(a['UNALIGNED_BASES']))+'</td><td>'+perc(a['UNALIGNED_BASES'],a['TOTAL_BASES'],1)+'</td></tr>'
  of.write(unaligned_bases_string+"\n")
  aligned_bases_string = '<tr><td>- Aligned bases</td><td>'+str(addcommas(a['ALIGNED_BASES']))+'</td><td>'+perc(a['ALIGNED_BASES'],a['TOTAL_BASES'],1)+'</td></tr>'
  of.write(aligned_bases_string+"\n")
  single_align_bases_string = '<tr><td>--- Single-aligned bases</td><td>'+str(addcommas(a['SINGLE_ALIGN_BASES']))+'</td><td>'+perc(a['SINGLE_ALIGN_BASES'],a['TOTAL_BASES'],1)+'</td></tr>'
  of.write(single_align_bases_string+"\n")
  gapped_align_bases_string = '<tr><td>--- Other-aligned bases</td><td>'+str(addcommas(a['GAPPED_ALIGN_BASES']))+'</td><td>'+perc(a['GAPPED_ALIGN_BASES'],a['TOTAL_BASES'],2)+'</td></tr>'
  of.write(gapped_align_bases_string+"\n")
  ostr = '''
    </table>
    <div id="legend_block">
      <table id="align_legend">
        <tr><td>Unaligned</td><td><div id="unaligned_leg"></div></td></tr>
        <tr><td>Trans-chimeric alignment</td><td><div id="chimeric_leg"></div></td></tr>
        <tr><td>Self-chimeric alignment</td><td><div id="selfchimeric_leg"></div></td></tr>
        <tr><td>Gapped alignment</td><td><div id="gapped_leg"></div></td></tr>
        <tr><td>Single alignment</td><td><div id="single_leg"></div></td></tr>
      </table>
    </div>
  </div>
  <div id="gapped_image_block">
    <div class="rhead">Summary [<a href="plots/alignments.pdf">pdf</a>]</div>
    <img id="gapped_image" src="plots/alignments.png">
  </div>   
</div>

<div class="clear"></div>
<hr>
<div class="subject_title">Coverage analysis &nbsp;&nbsp;&nbsp;&nbsp;<span class="highlight">'''
  of.write(ostr+"\n")
  of.write(perc(totcov,totref,2)+"\n")
  ostr = '''
  </span> <span class="highlight2">reference sequences covered</span>
</div>
<div id="coverage_block">
  <div class="rhead">Coverage of reference sequences [<a href="plots/covgraph.pdf">pdf<sub>1</sub></a>] [<a href="plots/perchrdepth.pdf">pdf<sub>2</sub></a>]</div>
  <div id="cov_images">
    <img id="covimg" class="square_image" src="plots/covgraph.png">
    <img id="chrcovimg" class="square_image" src="plots/perchrdepth.png">
  </div>
  <div class="clear"></div>
</div>
<hr>
<div class="subject_title"><table><tr><td class="c1">Rarefraction analysis</td><td class="c2"><span class="highlight">'''
  of.write(ostr)
  of.write(str(addcommas(geneany))+"\n")
  ostr = '''
  </span></td><td class="c3"><span class="highlight2">Genes detected</span></td><td class="c4"><span class="highlight">'''
  of.write(ostr)
  of.write(str(addcommas(genefull))+"\n")
  ostr = '''
  </span></td><td class="c5"><span class="highlight2">Full-length genes</span></td></tr></table>
</div>
<div id="annotated_block">
  <div class="insquare">
    <div class="rhead">Gene detection rarefraction [<a href="plots/gene_rarefraction.pdf">pdf</a>]</div>
    <img src="plots/gene_rarefraction.png">
  </div>
  <div class="insquare">
    <div class="rhead">Transcript detection rarefraction [<a href="plots/transcript_rarefraction.pdf">pdf</a>]</div>
    <img src="plots/transcript_rarefraction.png">
  </div>
  <div class="clear"></div>
  <div class="insquare">
    <table id="rarefraction_stats_table">
      <tr><td class="rhead" colspan="3">Rarefraction stats</td></tr>
      <tr class="bold"><td>Feature</td><td>Criteria</td><td>Count</td></tr>'''
  of.write(ostr+"\n")
  of.write('<tr><td>Gene</td><td>full-length</td><td>'+str(addcommas(genefull))+'</td></tr>')
  of.write('<tr><td>Gene</td><td>any match</td><td>'+str(addcommas(geneany))+'</td></tr>')
  of.write('<tr><td>Transcript</td><td>full-length</td><td>'+str(addcommas(txfull))+'</td></tr>')
  of.write('<tr><td>Transcript</td><td>any match</td><td>'+str(addcommas(txany))+'</td></tr>')
  of.write('<tr><td>Locus</td><td></td><td>'+str(addcommas(locuscount))+'</td></tr>')
  ostr='''
    </table>
    <table id="rarefraction_legend">
      <tr><td>Any match</td><td><div class="first"></div></td></tr>
      <tr><td>full-length</td><td><div class="second"></div></td></tr>
      <tr><td class="about" colspan="2">vertical line height indicates 5%-95% CI of simulation</td></tr>
    </table>
  </div>
  <div class="insquare">
    <div class="rhead">Locus detection rarefraction [<a href="plots/gene_rarefraction.pdf">pdf</a>]</div>
    <img id="locus_rarefraction_image" src="plots/locus_rarefraction.png">
  </div>
</div>
<div class="clear"></div>
<hr>
<div class="subject_title">Error pattern analysis &nbsp;&nbsp;&nbsp;&nbsp;<span class="highlight">'''
  of.write(ostr+"\n")
  error_rate = perc(e['ANY_ERROR'],e['ALIGNMENT_BASES'],3)
  of.write(error_rate)
  ostr='''
  </span> <span class="highlight2">error rate</span></div>
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
  ostr = '''
      <tr class="rhead"><td colspan="3">Base stats</td></tr>'''
  of.write(ostr+"\n")
  bases_analyzed_string = '<tr><td>Bases analyzed</td><td>'+str(addcommas(e['ALIGNMENT_BASES']))+'</td><td></td></tr>'
  of.write(bases_analyzed_string+"\n")
  correctly_aligned_string = '<tr><td>- Correctly aligned bases</td><td>'+str(addcommas(e['ALIGNMENT_BASES']-e['ANY_ERROR']))+'</td><td>'+perc((e['ALIGNMENT_BASES']-e['ANY_ERROR']),e['ALIGNMENT_BASES'],1)+'</td></tr>'
  of.write(correctly_aligned_string+"\n")
  total_error_string = '<tr><td>- Total error bases</td><td>'+str(addcommas(e['ANY_ERROR']))+'</td><td>'+perc(e['ANY_ERROR'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(total_error_string+"\n")
  mismatched_string = '<tr><td>--- Mismatched bases</td><td>'+str(addcommas(e['MISMATCHES']))+'</td><td>'+perc(e['MISMATCHES'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(mismatched_string+"\n")
  deletion_string = '<tr><td>--- Deletion bases</td><td>'+str(addcommas(e['ANY_DELETION']))+'</td><td>'+perc(e['ANY_DELETION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(deletion_string+"\n")
  complete_deletion_string = '<tr><td>----- Complete deletion bases</td><td>'+str(addcommas(e['COMPLETE_DELETION']))+'</td><td>'+perc(e['COMPLETE_DELETION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(complete_deletion_string+"\n")
  homopolymer_deletion_string = '<tr><td>----- Homopolymer deletion bases</td><td>'+str(addcommas(e['HOMOPOLYMER_DELETION']))+'</td><td>'+perc(e['HOMOPOLYMER_DELETION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(homopolymer_deletion_string+"\n")
  insertion_string = '<tr><td>--- Insertion bases</td><td>'+str(addcommas(e['ANY_INSERTION']))+'</td><td>'+perc(e['ANY_INSERTION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(insertion_string+"\n")
  complete_insertion_string = '<tr><td>----- Complete insertion bases</td><td>'+str(addcommas(e['COMPLETE_INSERTION']))+'</td><td>'+perc(e['COMPLETE_INSERTION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(complete_insertion_string+"\n")
  homopolymer_insertion_string = '<tr><td>----- Homopolymer insertion bases</td><td>'+str(addcommas(e['HOMOPOLYMER_INSERTION']))+'</td><td>'+perc(e['HOMOPOLYMER_INSERTION'],e['ALIGNMENT_BASES'],3)+'</td></tr>'
  of.write(homopolymer_insertion_string+"\n")
  ostr = '''
    </table>
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

def addcommas(val):
  return locale.format("%d",val,grouping=True)

if __name__=="__main__":
  main()
