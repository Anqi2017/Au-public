#!/usr/bin/python
import argparse, sys, os, time, re, gzip, locale, inspect
from subprocess import Popen, PIPE

# BAM imports
import bam_traversal
import get_platform_report
import gpd_loci_analysis
import gpd_to_exon_distro
import make_alignment_plot
import depth_to_coverage_report
import locus_bed_to_rarefraction

# BAM + reference imports
import bam_to_context_error_plot
import bam_to_alignment_error_plot

# BAM + annotation
import  annotate_from_genomic_features
import  get_depth_subset
import  annotated_length_analysis
import  gpd_annotation_to_rarefraction
import  annotated_read_bias_analysis

#bring in the folder to the path for our utilities
pythonfolder_loc = "../../../utilities"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)

# BAM
import gpd_to_bed_depth
import genepred_to_bed

# BAM + annotation
import  gpd_annotate

# read count
rcnt = -1

def main(args):

  if not os.path.exists(args.tempdir+'/plots'):
    os.makedirs(args.tempdir+'/plots')
  if not os.path.exists(args.tempdir+'/data'):
    os.makedirs(args.tempdir+'/data')
  if not os.path.exists(args.tempdir+'/logs'):
    os.makedirs(args.tempdir+'/logs')

  ## Extract data that can be realized from the bam
  make_data_bam(args)

  ## Extract data that can be realized from the bam and reference
  if args.reference:
    make_data_bam_reference(args)

  ## Extract data that can be realized from bam and reference annotation
  if args.annotation:
    make_data_bam_annotation(args)

  # Write params file
  of = open(args.tempdir+'/data/params.txt','w')
  for arg in vars(args):
    of.write(arg+"\t"+str(getattr(args,arg))+"\n")
  of.close()

def make_data_bam(args):
  # Get the data necessary for making tables and reports

  # 1. Traverse bam file to describe alignment details
  udir = os.path.dirname(os.path.realpath(__file__))
  cmd = udir+'/bam_traversal.py '+args.input+' -o '+args.tempdir+'/data/ '
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
  # Calling the command through python
  bam_traversal.external_cmd(cmd)

  # Now we can find any known reads
  # 2. Go through read names to find if there are platform-specific read names present
  sys.stderr.write("Can we find any known read types\n")
  cmd = udir+'/get_platform_report.py '+args.tempdir+'/data/lengths.txt.gz '
  cmd += args.tempdir+'/data/special_report'
  sys.stderr.write(cmd+"\n") 
  get_platform_report.external_cmd(cmd)

  # Check for pacbio to see if we need to make a graph for it
  do_pb = False
  with open(args.tempdir+'/data/special_report') as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      if f[0]=='PB': 
        do_pb = True
        break
  if do_pb:
    cmd = 'Rscript '+udir+'/plot_pacbio.r '+args.tempdir+'/data/special_report.pacbio '+args.tempdir+'/plots/pacbio.png'
    sys.stderr.write(cmd+"\n")
    mycall(cmd,args.tempdir+'/logs/special_report_pacbio_png')
    cmd = 'Rscript '+udir+'/plot_pacbio.r '+args.tempdir+'/data/special_report.pacbio '+args.tempdir+'/plots/pacbio.pdf'
    sys.stderr.write(cmd+"\n")
    mycall(cmd,args.tempdir+'/logs/special_report_pacbio_pdf')

  # 3. Go through the genepred file and get a depth bed for our best alignments
  sys.stderr.write("Go through genepred best alignments and make a bed depth file\n")
  cmd = "gpd_to_bed_depth.py "+args.tempdir+'/data/best.sorted.gpd.gz -o '+args.tempdir+'/data/depth.sorted.bed.gz'
  sys.stderr.write("Generate the depth bed for the mapped reads\n")
  sys.stderr.write(cmd+"\n")
  gpd_to_bed_depth.external_cmd(cmd)

  global rcnt #read count
  rcnt = 0
  tinf = gzip.open(args.tempdir+'/data/lengths.txt.gz')
  for line in tinf:  rcnt += 1
  tinf.close()

  # For now reporting loci will be optional until it can be tested and optimized.
  if args.do_loci:
    # 4. Go through the best alignments and look for loci
    sys.stderr.write("Approximate loci and mapped read distributions among them.\n")
    cmd = udir+"/gpd_loci_analysis.py "+args.tempdir+'/data/best.sorted.gpd.gz -o '+args.tempdir+'/data/loci-all.bed.gz --output_loci '+args.tempdir+'/data/loci.bed.gz'
    cmd += ' --downsample '+str(args.locus_downsample)+' '
    cmd += ' --threads '+str(args.threads)+' '
    if args.min_depth:
      cmd += ' --min_depth '+str(args.min_depth)
    if args.min_depth:
      cmd += ' --min_coverage_at_depth '+str(args.min_coverage_at_depth)
    if args.min_exon_count:
      cmd += ' --min_exon_count '+str(args.min_exon_count)
    sys.stderr.write(cmd+"\n")
    gpd_loci_analysis.external_cmd(cmd)

    cmd = udir+"/locus_bed_to_rarefraction.py "+args.tempdir+'/data/loci.bed.gz -o '+args.tempdir+'/data/locus_rarefraction.txt'
    cmd += ' --threads '+str(args.threads)+' '
    cmd += ' --original_read_count '+str(rcnt)+' '
    sys.stderr.write("Make rarefraction curve\n")
    sys.stderr.write(cmd+"\n")
    locus_bed_to_rarefraction.external_cmd(cmd)

    sys.stderr.write("Make locus rarefraction plot\n")
    for ext in ['png','pdf']:
      cmd = 'Rscript '+udir+'/plot_annotation_rarefractions.r '+\
               args.tempdir+'/plots/locus_rarefraction.'+ext+' '+\
               'locus'+' '+\
               args.tempdir+'/data/locus_rarefraction.txt '+\
               '#FF000088 '
      sys.stderr.write(cmd+"\n")
      mycall(cmd,args.tempdir+'/logs/plot_locus_rarefraction_'+ext)

  # 6. Alignment plot preparation
  sys.stderr.write("Get ready for alignment plot\n")
  cmd = udir+'/make_alignment_plot.py '+args.tempdir+'/data/lengths.txt.gz '
  cmd += ' --output_stats '+args.tempdir+'/data/alignment_stats.txt '
  cmd += ' --output '+args.tempdir+'/plots/alignments.png '
  cmd += args.tempdir+'/plots/alignments.pdf'
  sys.stderr.write("Make alignment plots\n")
  sys.stderr.write(cmd+"\n")
  make_alignment_plot.external_cmd(cmd)

  # 7. Make depth reports
  sys.stderr.write("Making depth reports\n")
  cmd = udir+'/depth_to_coverage_report.py '+args.tempdir+'/data/depth.sorted.bed.gz '+args.tempdir+'/data/chrlens.txt -o '+args.tempdir+'/data'
  sys.stderr.write(cmd+"\n")
  depth_to_coverage_report.external_cmd(cmd)

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

  #Get the exon distribution
  sys.stderr.write("Get the exon distributions\n")
  cmd = udir+'/gpd_to_exon_distro.py '
  cmd += args.tempdir+'/data/best.sorted.gpd.gz -o '
  cmd += args.tempdir+'/data/exon_size_distro.txt.gz'
  sys.stderr.write(cmd+"\n")
  gpd_to_exon_distro.external_cmd(cmd)

  cmd = 'Rscript '+udir+'/plot_exon_distro.r '+args.tempdir+'/data/exon_size_distro.txt.gz '+args.tempdir+'/plots/exon_size_distro.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/exon_size_distro_png')
  cmd = 'Rscript '+udir+'/plot_exon_distro.r '+args.tempdir+'/data/exon_size_distro.txt.gz '+args.tempdir+'/plots/exon_size_distro.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/exon_size_distro_pdf')

  # Make a UCSC compatible bed file
  sys.stderr.write("Make a UCSC genome browser compatible bed file\n")
  cmd = 'genepred_to_bed.py --headername '+args.input+':best '
  cmd += ' '+args.tempdir+'/data/best.sorted.gpd.gz'
  cmd += ' -o '+args.tempdir+'/data/best.sorted.bed.gz'
  cmd += ' --color red'
  sys.stderr.write(cmd+"\n")
  genepred_to_bed.external_cmd(cmd)

  cmd = 'genepred_to_bed.py --headername '+args.input+':trans-chimera '
  cmd += ' '+args.tempdir+'/data/chimera.gpd.gz'
  cmd += ' -o '+args.tempdir+'/data/chimera.bed.gz'
  cmd += ' --color blue'
  sys.stderr.write(cmd+"\n")
  genepred_to_bed.external_cmd(cmd)

  cmd = 'genepred_to_bed.py --headername '+args.input+':gapped '
  cmd += ' '+args.tempdir+'/data/gapped.gpd.gz'
  cmd += ' -o '+args.tempdir+'/data/gapped.bed.gz'
  cmd += ' --color orange'
  sys.stderr.write(cmd+"\n")
  genepred_to_bed.external_cmd(cmd)

  cmd = 'genepred_to_bed.py --headername '+args.input+':self-chimera '
  cmd += ' '+args.tempdir+'/data/technical_chimeras.gpd.gz'
  cmd += ' -o '+args.tempdir+'/data/technical_chimeras.bed.gz'
  cmd += ' --color green'
  sys.stderr.write(cmd+"\n")
  genepred_to_bed.external_cmd(cmd)

  cmd = 'genepred_to_bed.py --headername '+args.input+':self-atypical '
  cmd += ' '+args.tempdir+'/data/technical_atypical_chimeras.gpd.gz'
  cmd += ' -o '+args.tempdir+'/data/technical_atypical_chimeras.bed.gz'
  cmd += ' --color purple'
  sys.stderr.write(cmd+"\n")
  genepred_to_bed.external_cmd(cmd)

def make_data_bam_reference(args):
  # make the context error plots
  udir = os.path.dirname(os.path.realpath(__file__))

  # 1. Context error
  cmd = udir+'/bam_to_context_error_plot.py '+args.input+' -r '+args.reference+' --target --output_raw '+args.tempdir+'/data/context_error_data.txt -o '+args.tempdir+'/plots/context_plot.png '+args.tempdir+'/plots/context_plot.pdf'
  if args.context_error_scale:
    cmd += ' --scale '+' '.join([str(x) for x in args.context_error_scale])
  if args.context_error_stopping_point:
    cmd += ' --stopping_point '+str(args.context_error_stopping_point)
  sys.stderr.write("Making context plot\n")
  sys.stderr.write(cmd+"\n")
  bam_to_context_error_plot.external_cmd(cmd)

  # 2. Alignment overall error
  cmd = udir+'/bam_to_alignment_error_plot.py '+args.input+' -r '+args.reference+' --output_stats '+args.tempdir+'/data/error_stats.txt --output_raw '+args.tempdir+'/data/error_data.txt -o '+args.tempdir+'/plots/alignment_error_plot.png '+args.tempdir+'/plots/alignment_error_plot.pdf'
  if args.alignment_error_scale:
    cmd += ' --scale '+' '.join([str(x) for x in args.alignment_error_scale])
  if args.alignment_error_max_length:
    cmd += ' --max_length '+str(args.alignment_error_max_length)
  sys.stderr.write("Making alignment error plot\n")
  sys.stderr.write(cmd+"\n")
  bam_to_alignment_error_plot.external_cmd(cmd)

def make_data_bam_annotation(args):
  udir = os.path.dirname(os.path.realpath(__file__))

  # 1. Use annotations to identify genomic features (Exon, Intron, Intergenic)
  # And assign membership to reads
  # Stores the feature bed files in a beds folder
  cmd = udir+'/annotate_from_genomic_features.py --output_beds '+args.tempdir+'/data/beds '
  cmd += args.tempdir+'/data/best.sorted.gpd.gz '+args.annotation+' '
  cmd += args.tempdir+'/data/chrlens.txt -o '+args.tempdir+'/data/read_genomic_features.txt.gz'
  sys.stderr.write("Finding genomic features and assigning reads membership\n")
  sys.stderr.write(cmd+"\n")
  annotate_from_genomic_features.external_cmd(cmd)

  # 2. Get depth distributions for each feature subset
  # now get depth subsets
  sys.stderr.write("get depths of features\n")
  cmd = udir+'/get_depth_subset.py '+args.tempdir+'/data/depth.sorted.bed.gz '
  cmd += args.tempdir+'/data/beds/exon.bed -o '
  cmd += args.tempdir+'/data/exondepth.bed.gz'
  sys.stderr.write(cmd+"\n")
  get_depth_subset.external_cmd(cmd)
  
  cmd = udir+'/get_depth_subset.py '+args.tempdir+'/data/depth.sorted.bed.gz '
  cmd += args.tempdir+'/data/beds/intron.bed -o '
  cmd += args.tempdir+'/data/introndepth.bed.gz'
  sys.stderr.write(cmd+"\n")
  get_depth_subset.external_cmd(cmd)

  cmd = udir+'/get_depth_subset.py '+args.tempdir+'/data/depth.sorted.bed.gz '
  cmd += args.tempdir+'/data/beds/intergenic.bed -o '
  cmd += args.tempdir+'/data/intergenicdepth.bed.gz'
  sys.stderr.write(cmd+"\n")
  get_depth_subset.external_cmd(cmd)

  # 3. Plot the feature depth
  cmd = 'Rscript '+udir+'/plot_feature_depth.r '
  cmd += args.tempdir+'/data/depth.sorted.bed.gz '
  cmd += args.tempdir+'/data/exondepth.bed.gz '
  cmd += args.tempdir+'/data/introndepth.bed.gz '
  cmd += args.tempdir+'/data/intergenicdepth.bed.gz '
  cmd += args.tempdir+'/plots/feature_depth.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/featuredepth_png')  

  cmd = 'Rscript '+udir+'/plot_feature_depth.r '
  cmd += args.tempdir+'/data/depth.sorted.bed.gz '
  cmd += args.tempdir+'/data/exondepth.bed.gz '
  cmd += args.tempdir+'/data/introndepth.bed.gz '
  cmd += args.tempdir+'/data/intergenicdepth.bed.gz '
  cmd += args.tempdir+'/plots/feature_depth.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/featuredepth_pdf')  

  # 4. Generate plots from reads assigend to features
  sys.stderr.write("Plot read assignment to genomic features\n")
  cmd = 'Rscript '+udir+'/plot_annotated_features.r '
  cmd += args.tempdir+'/data/read_genomic_features.txt.gz '
  cmd += args.tempdir+'/plots/read_genomic_features.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/read_genomic_features_png')  

  cmd = 'Rscript '+udir+'/plot_annotated_features.r '
  cmd += args.tempdir+'/data/read_genomic_features.txt.gz '
  cmd += args.tempdir+'/plots/read_genomic_features.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/read_genomic_features_pdf')  
  
  # 5. annotated the best mappend read mappings
  cmd = 'gpd_annotate.py '+args.tempdir+'/data/best.sorted.gpd.gz -r '+args.annotation+' -o '+args.tempdir+'/data/annotbest.txt.gz'
  if args.threads:
    cmd += ' --threads '+str(args.threads)
  sys.stderr.write("Annotating reads\n")
  sys.stderr.write(cmd+"\n")
  gpd_annotate.external_cmd(cmd)

  # 6. Make plots of the transcript lengths
  sys.stderr.write("Make plots from transcript lengths\n")
  cmd = 'Rscript '+udir+'/plot_transcript_lengths.r '
  cmd += args.tempdir+'/data/annotbest.txt.gz '
  cmd += args.tempdir+'/plots/transcript_distro.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/transcript_distro_png')  
  
  sys.stderr.write("Make plots from transcript lengths\n")
  cmd = 'Rscript '+udir+'/plot_transcript_lengths.r '
  cmd += args.tempdir+'/data/annotbest.txt.gz '
  cmd += args.tempdir+'/plots/transcript_distro.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/transcript_distro_pdf')  
  

  # 7. Make length distributions for plotting
  sys.stderr.write("making length distributions from annotations\n")
  cmd = udir+'/annotated_length_analysis.py '
  cmd += args.tempdir+'/data/best.sorted.gpd.gz '
  cmd += args.tempdir+'/data/annotbest.txt.gz '
  cmd += '-o '+args.tempdir+'/data/annot_lengths.txt.gz'
  sys.stderr.write(cmd+"\n")
  annotated_length_analysis.external_cmd(cmd)

  # 8. Plot length distributions
  cmd = 'Rscript '+udir+'/plot_annotation_analysis.r '
  cmd += args.tempdir+'/data/annot_lengths.txt.gz '
  cmd += args.tempdir+'/plots/annot_lengths.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/annot_lengths_png')  

  cmd = 'Rscript '+udir+'/plot_annotation_analysis.r '
  cmd += args.tempdir+'/data/annot_lengths.txt.gz '
  cmd += args.tempdir+'/plots/annot_lengths.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/annot_lengths_pdf')  

  # 9. Get rarefraction curve data
  global rcnt
  if rcnt < 0:
    sys.stderr.write("Getting read count\n")
    rcnt = 0
    tinf = gzip.open(args.tempdir+'/data/lengths.txt.gz')
    for line in tinf:  rcnt += 1
    tinf.close()
  sys.stderr.write("Writing rarefraction curves\n")
  cmd =  udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --samples_per_xval '+str(args.samples_per_xval)
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --gene -o '+args.tempdir+'/data/gene_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)

  cmd =  udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --samples_per_xval '+str(args.samples_per_xval)
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --transcript -o '+args.tempdir+'/data/transcript_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)

  cmd =  udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --samples_per_xval '+str(args.samples_per_xval)
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --full --gene -o '+args.tempdir+'/data/gene_full_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)

  cmd =  udir+'/gpd_annotation_to_rarefraction.py '+args.tempdir+'/data/annotbest.txt.gz '
  cmd += ' --samples_per_xval '+str(args.samples_per_xval)
  cmd += ' --original_read_count '+str(rcnt)
  cmd += ' --threads '+str(args.threads)
  cmd += ' --full --transcript -o '+args.tempdir+'/data/transcript_full_rarefraction.txt'
  sys.stderr.write(cmd+"\n")
  gpd_annotation_to_rarefraction.external_cmd(cmd)

  # 10. Plot the rarefraction curves
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

  # 11. Use annotation outputs to check for  bias
  sys.stderr.write("Prepare bias data\n")
  cmd = udir+'/annotated_read_bias_analysis.py '+\
        args.tempdir+'/data/best.sorted.gpd.gz '+\
        args.annotation+' '+ args.tempdir+'/data/annotbest.txt.gz '+\
        '-o '+args.tempdir+'/data/bias_table.txt.gz '+\
        '--output_counts '+args.tempdir+'/data/bias_counts.txt'
  sys.stderr.write(cmd+"\n")
  annotated_read_bias_analysis.external_cmd(cmd)

  # 12. Plot bias
  cmd = 'Rscript '+udir+'/plot_bias.r '+args.tempdir+'/data/bias_table.txt.gz '+\
        args.tempdir+'/plots/bias.png'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/bias_png.log')

  cmd = 'Rscript '+udir+'/plot_bias.r '+args.tempdir+'/data/bias_table.txt.gz '+\
        args.tempdir+'/plots/bias.pdf'
  sys.stderr.write(cmd+"\n")
  mycall(cmd,args.tempdir+'/logs/bias_pdf.log')

  return


def mycall(cmd,lfile):
  ofe = open(lfile+'.err','w')
  ofo = open(lfile+'.out','w')
  p = Popen(cmd.split(),stderr=ofe,stdout=ofo)
  p.communicate()
  ofe.close()
  ofo.close()
  return

#def do_inputs():
#  # Setup command line inputs
#  parser=argparse.ArgumentParser(description="Create an output report",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#  parser.add_argument('input',help="INPUT FILE or '-' for STDIN")
#  parser.add_argument('-o','--output',help="OUTPUT Folder or STDOUT if not set")
#  parser.add_argument('--portable_output',help="OUTPUT file in a portable html format")
#  group1 = parser.add_mutually_exclusive_group(required=True)
#  group1.add_argument('-r','--reference',help="Reference Fasta")
#  group1.add_argument('--no_reference',action='store_true',help="No Reference Fasta")
#  parser.add_argument('--annotation',help="Reference annotation genePred")
#  parser.add_argument('--threads',type=int,default=1,help="INT number of threads to run. Default is system cpu count")
#  # Temporary working directory step 1 of 3 - Definition
#  parser.add_argument('--tempdir',required=True,help="This temporary directory will be used, but will remain after executing.")
#
#  ### Parameters for alignment plots
#  parser.add_argument('--min_aligned_bases',type=int,default=50,help="for analysizing alignment, minimum bases to consider")
#  parser.add_argument('--max_query_overlap',type=int,default=10,help="for testing gapped alignment advantage")
#  parser.add_argument('--max_target_overlap',type=int,default=10,help="for testing gapped alignment advantage")
#  parser.add_argument('--max_query_gap',type=int,help="for testing gapped alignment advantge")
#  parser.add_argument('--max_target_gap',type=int,default=500000,help="for testing gapped alignment advantage")
#  parser.add_argument('--required_fractional_improvement',type=float,default=0.2,help="require gapped alignment to be this much better (in alignment length) than single alignment to consider it.")
#  
#  ### Parameters for locus analysis
#  parser.add_argument('--do_loci',action='store_true',help="This analysis is time consuming at the moment so don't do it unless necessary")
#  parser.add_argument('--min_depth',type=float,default=1.5,help="require this or more read depth to consider locus")
#  parser.add_argument('--min_coverage_at_depth',type=float,default=0.8,help="require at leas this much of the read be covered at min_depth")
#  parser.add_argument('--min_exon_count',type=int,default=2,help="Require at least this many exons in a read to consider assignment to a locus")
#  parser.add_argument('--locus_downsample',type=int,default=100,help="Only include up to this many long reads in a locus\n")
#
#  ### Params for alignment error plot
#  parser.add_argument('--alignment_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
#  parser.add_argument('--alignment_error_max_length',type=int,default=100000,help="The maximum number of alignment bases to calculate error from")
#  
#  ### Params for context error plot
#  parser.add_argument('--context_error_scale',nargs=6,type=float,help="<ins_min> <ins_max> <mismatch_min> <mismatch_max> <del_min> <del_max>")
#  parser.add_argument('--context_error_stopping_point',type=int,default=1000,help="Sample at least this number of each context")
#
#  ## Params for rarefraction plots
#  parser.add_argument('--samples_per_xval',type=int,default=500)
#
#  args = parser.parse_args()
#  # Temporary working directory step 2 of 3 - Creation
#  setup_tempdir(args)
#  return args

def setup_tempdir(args):
  if not os.path.exists(args.tempdir):
    os.makedirs(args.tempdir.rstrip('/'))
  if not os.path.exists(args.tempdir):
    sys.stderr.write("ERROR: Problem creating temporary directory\n")
    sys.exit()
  return 

def external(args):
  main(args)

if __name__=="__main__":
  sys.stderr.write("excute as prepare all data as main\n")
  #do our inputs
  # Can disable calling as main
  #args = do_inputs()
  #main(args)
