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
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default='/tmp',help='FOLDERNAME location of random temporary directory')
  group.add_argument('--specific_tempdir',help='FOLDERNAME location of actual temporary directory. will not remove during cleanup.')
  parser.add_argument('--output',default='output_pre-IDP_step-1_from_raw',help='FOLDERNAME of output must not already exist')
  parser.add_argument('--ccs_hq_acc',type=int,default=95,help='INT accuracy of high quality ccs reads')
  parser.add_argument('--ccs_hq_passes',type=int,default=2,help='INT number of passes for high quality ccs reads')
  parser.add_argument('--ccs_lq_acc',type=int,default=75,help='INT accuracy of low quality ccs reads')
  parser.add_argument('--ccs_lq_passes',type=int,default=0,help='INT number of passes for low quality ccs reads')
  parser.add_argument('--subreads_acc',type=int,default=75,help='INT minimum accuracy of subreads')
  parser.add_argument('--save_tempdir',help='DIRECTORYNAME name of a directory to be created that has the temporary folder in it.  Unnecessary if you already use a specific_tempdir')
  args = parser.parse_args()

  if args.threads==0:
    args.threads = multiprocessing.cpu_count()

  if not os.path.isdir(args.smrtanalysis_path.rstrip('/')+'/current'):
    sys.stderr.write("ERROR could not locate folder current in smrtanalysis_path: "+args.smrtanalysis_path+"\n")
    sys.exit()
  if os.path.isdir(args.output):
    sys.stderr.write("ERROR "+args.output+" output folder must not already exist.\n")
    sys.exit()
  if args.save_tempdir:
    if os.path.isdir(args.save_tempdir):
      sys.stderr.write("ERROR "+args.save_tempdir+" folder must not already exist.\n")
      sys.exit()
  args.smrtanalysis_path = args.smrtanalysis_path.rstrip('/')

  tdir = setup_temporary_directory(args)  
  sys.stderr.write("working in "+tdir+"\n")

  if not os.path.exists(tdir+'/output'):
    os.makedirs(tdir+'/output')

  of_log = open(tdir+'/output/LOG','w')
  of_log.write(str(sys.argv)+"\n")

  sys.stderr.write("Extracting high quality ccs reads from "+args.pacbio_raw+"\n")
  execute_ccs(tdir,args,args.ccs_hq_passes,args.ccs_hq_acc,'ccs95_out')
  of_log.write("High quality ccs read accuracy: "+str(args.ccs_hq_acc)+"\n")
  of_log.write("High quality ccs read passes: "+str(args.ccs_hq_passes)+"\n")
  [ccs_hq_fasta,ccs_hq_fastq,ccs_hq_h5] = get_ccs_paths(tdir,'ccs95_out')
  if not os.path.exists(tdir+'/output/ccs_hq'):
    os.makedirs(tdir+'/output/ccs_hq')
  copyfile(ccs_hq_fasta,tdir+'/output/ccs_hq/ccs_hq.fa')
  copyfile(ccs_hq_fastq,tdir+'/output/ccs_hq/ccs_hq.fq')
  copyfile(ccs_hq_h5,tdir+'/output/ccs_hq/ccs_hq.h5')
  
  sys.stderr.write("Extracting low quality ccs reads from "+args.pacbio_raw+"\n")
  execute_ccs(tdir,args,args.ccs_lq_passes,args.ccs_lq_acc,'ccs75_out')
  of_log.write("Low quality ccs read accuracy: "+str(args.ccs_lq_acc)+"\n")
  of_log.write("Low quality ccs read passes: "+str(args.ccs_lq_passes)+"\n")
  [ccs_lq_fasta,ccs_lq_fastq,ccs_lq_h5] = get_ccs_paths(tdir,'ccs75_out')
  if not os.path.exists(tdir+'/output/ccs_lq'):
    os.makedirs(tdir+'/output/ccs_lq')
  copyfile(ccs_lq_fasta,tdir+'/output/ccs_lq/ccs_lq.fa')
  copyfile(ccs_lq_fastq,tdir+'/output/ccs_lq/ccs_lq.fq')
  copyfile(ccs_lq_h5,tdir+'/output/ccs_lq/ccs_lq.h5')

  sys.stderr.write("Extracting subreads from "+args.pacbio_raw+"\n")
  execute_subreads(tdir,args,args.subreads_acc,'subreads_out')
  of_log.write("Subread read accuracy: "+str(args.subreads_acc)+"\n")
  [subreads_fasta,subreads_fastq] = get_subreads_paths(tdir,'subreads_out')
  if not os.path.exists(tdir+'/output/subreads'):
    os.makedirs(tdir+'/output/subreads')
  copyfile(ccs_lq_fasta,tdir+'/output/subreads/subreads.fa')
  copyfile(ccs_lq_fastq,tdir+'/output/subreads/subreads.fq')

  sys.stderr.write("Get a set of reads to correct 75-95 and longest subreads\n")
  to_correct_fasta = tdir+'/output/ccs_lq_and_longest_subreads_to_correct.fa'
  get_sequences_to_correct(ccs_hq_fasta,ccs_lq_fasta,subreads_fasta,to_correct_fasta)

  # make a set of non redundant uncorrected reads
  of_un = open(tdir+'/output/lr_nonredundant_uncorrected.fa','w')
  with open(ccs_hq_fasta) as inf:
    for line in inf:
      of_un.write(line)
  with open(tdir+'/output/ccs_lq_and_longest_subreads_to_correct.fa') as inf:
    for line in inf:
      of_un.write(line)
  of_un.close()

  of_log.close()
  if args.save_tempdir:
    copytree(tdir,args.save_tempdir)
  copytree(tdir+'/output',args.output)
  if not args.specific_tempdir:
    rmtree(tdir)

def get_sequences_to_correct(ccs_hq_file,ccs_lq_file,subreads_file,output_fasta):
  ccs95 = read_fasta_into_array(ccs_hq_file)
  #get the names and base names from ccs95
  ccs95basenames = set()
  prog = re.compile('(^[^\/]+\/\d+)\/')
  of = open(output_fasta,'w')
  for entry in ccs95:
    m = prog.match(entry['name'])
    if not m: 
      sys.stderr.write('ERROR trouble parsing name '+entry['name']+"\n")
      sys.exit()
    basename = m.group(1)
    ccs95basenames.add(basename)
  ccs95 = []
  ccs90 = read_fasta_into_array(ccs_lq_file)
  ccs90to95basenames = set()
  for entry in ccs90:
    m = prog.match(entry['name'])
    if not m:
      sys.stderr.write('trouble parsing name'+entry['name']+"\n")
      sys.exit()
    basename = m.group(1)
    if basename in ccs95basenames: continue
    of.write('>'+entry['name']+"\n")
    of.write(entry['seq']+"\n")
    ccs90to95basenames.add(basename)
  ccs90 = []
  sub75 = read_fasta_into_array(subreads_file)
  sub75lengths = {}
  for entry in sub75:
    m = prog.match(entry['name'])
    if not m:
      sys.stderr.write('trouble parsing name '+entry['name']+"\n")
      sys.exit()
    basename = m.group(1)
    if basename in ccs95basenames: continue
    if basename in ccs90to95basenames: continue
    if basename not in sub75lengths: 
      sub75lengths[basename] = len(entry['seq'])
  printsub75 = {}
  for entry in sub75:
    m = prog.match(entry['name'])
    if not m:
      sys.stderr.write('trouble parsing name '+entry['name']+"\n")
      sys.exit()
    basename = m.group(1)
    if basename in ccs95basenames: continue
    if basename in ccs90to95basenames: continue
    if len(entry['seq']) == sub75lengths[basename]:
      printsub75[basename] = entry
  for basename in printsub75:
    entry = printsub75[basename]
    of.write('>'+entry['name']+"\n")
    of.write(entry['seq']+"\n")
  of.close()

def get_subreads_paths(tdir,output_dir):
  fasta = None
  fastq = None
  for filename in os.listdir(tdir+'/'+output_dir.rstrip('/')+'/data'):
    if re.search('\.fasta$',filename): 
      fasta = tdir+'/'+output_dir.rstrip('/')+'/data/'+filename
    elif re.search('\.fastq$',filename):
      fastq = tdir+'/'+output_dir.rstrip('/')+'/data/'+filename
  return [fasta, fastq]

def get_ccs_paths(tdir,output_dir):
  fasta = None
  fastq = None
  h5 = None
  for filename in os.listdir(tdir+'/'+output_dir):
    if re.search('\.fasta$',filename): 
      fasta = tdir+'/'+output_dir.rstrip('/')+'/'+filename
    elif re.search('\.fastq$',filename):
      fastq = tdir+'/'+output_dir.rstrip('/')+'/'+filename
    elif re.search('\.h5$',filename):
      h5 = tdir+'/'+output_dir.rstrip('/')+'/'+filename
  return [fasta, fastq, h5]

def execute_subreads(tdir,args,min_acc,output_dir):
  input_xml = '''<?xml version="1.0"?>
<pacbioAnalysisInputs>
  <dataReferences>
    <url ref="run:0000000-0000"><location>'''+os.path.abspath(args.pacbio_raw)+'''</location></url>
  </dataReferences>
</pacbioAnalysisInputs>'''
  of = open(tdir+'/input.xml','w')
  of.write(input_xml+"\n")
  of.close()

  settings_xml = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<smrtpipeSettings>
    <protocol version="2.3.0" id="RS_Subreads.1" editable="false">
    <application>Data Prep</application>
        <param name="name" label="Protocol Name">
            <value>RS_Subreads</value>
            <input type="text"/>
            <rule required="true" message="Protocol name is required"/>
        </param>
        <param name="description">
            <value>Filter subreads based on read length and quality, optionally splitting by barcode. Output FASTA and bas.h5 file of subreads.</value>
            <textarea></textarea>
        </param>
        <param name="version" hidden="true">
            <value>1</value>
            <input type="text"/>
            <rule type="digits" required="true" min="1.0"/>
        </param>
        <param name="state">
            <value>active</value>
            <input value="active" type="radio"/>
            <input value="inactive" type="radio"/>
        </param>
        <param name="control" hidden="true">
            <value></value>
        </param>
        <param name="fetch" hidden="true">
            <value>common/protocols/preprocessing/Fetch.1.xml</value>
        </param>
        <param name="filtering">
            <value>common/protocols/filtering/SFilter.1.xml</value>
            <select multiple="true">
                <import extension="xml" contentType="text/directory">common/protocols/filtering</import>
            </select>
        </param>
        <param name="barcode" editableInJob="true">
            <value>common/protocols/barcode/NoBarcode.1.xml</value>
            <select multiple="false">
                <import extension="xml" contentType="text/directory">common/protocols/barcode</import>
            </select>
        </param>
    </protocol>
    <moduleStage name="fetch" editable="true">
        <module label="Fetch v1" id="P_Fetch" editableInJob="true">
            <description>Sets up inputs</description>
        </module>
    </moduleStage>
    <moduleStage name="filtering" editable="true">
        <module label="SFilter v1" id="P_Filter" editableInJob="true">
            <description>This module filters reads based on a minimum subread length, polymerase read quality and polymerase read length.</description>
            <param name="minSubReadLength" label="Minimum Subread Length">
                <value>50</value>
                <title>Subreads shorter than this value (in base pairs) are filtered out and excluded from analysis.</title>
                <input type="text" size="3"/>
                <rule type="number" min="0.0" message="Value must be a positive integer"/>
            </param>
            <param name="readScore" label="Minimum Polymerase Read Quality">
                <value>'''+str(min_acc)+'''</value>
                <title>Polymerase reads with lower quality than this value are filtered out and excluded from analysis.</title>
                <input type="text"/>
                <rule type="number" min="0.0" message="Value must be between 0 and 100" max="100.0"/>
            </param>
            <param name="minLength" label="Minimum Polymerase Read Length">
                <value>50</value>
                <title>Polymerase reads shorter than this value (in base pairs) are filtered out and excluded from analysis.</title>
                <input type="text" size="3"/>
                <rule type="number" min="0.0" message="Value must be a positive integer"/>
            </param>
        </module>
        <module label="SFilter Reports v1" id="P_FilterReports" editableInJob="false"/>
    </moduleStage>
    <moduleStage name="barcode" editable="true"/>
    <fileName>RS_Subreads.1.xml</fileName>
</smrtpipeSettings>'''
  of = open(tdir+'/settings.xml','w')
  of.write(settings_xml+"\n")
  of.close()
  if not os.path.exists(tdir+'/mytemp'):
    os.makedirs(tdir+'/mytemp')
  if not os.path.exists(tdir+'/'+output_dir):
    os.makedirs(tdir+'/'+output_dir)
  cmd1  = '. '+args.smrtanalysis_path+'/current/etc/setup.sh && '
  cmd1 += args.smrtanalysis_path+'/current/analysis/bin/smrtpipe.py '
  cmd1 += '-D TMP='+tdir+'/mytemp -D SHARED_DIR='+tdir+'/mytemp '
  cmd1 += '--output='+tdir+'/'+output_dir+' '
  cmd1 += '--params='+tdir+'/settings.xml xml:'+tdir+'/input.xml'
  subprocess.call(cmd1,shell=True)
  return 

def execute_ccs(tdir,args,min_pass,min_acc,output_dir):
  of = open(tdir+'/rawname.txt','w')
  of.write(os.path.abspath(args.pacbio_raw)+"\n")
  of.close()
  setup_sh_file = args.smrtanalysis_path.rstrip('/')+'/current/etc/setup.sh'
  consensus_tools_sh_file = args.smrtanalysis_path.rstrip('/')+'/current/analysis/bin/ConsensusTools.sh'
  parameter_file = args.smrtanalysis_path.rstrip('/')+'/current/analysis/etc/algorithm_parameters/2014-09 '
  cmd1  = '. '+setup_sh_file+' && '+consensus_tools_sh_file+' CircularConsensus '
  cmd1 += '--minFullPasses '+str(min_pass)+' --minPredictedAccuracy '+str(min_acc)+' '
  cmd1 += '--parameters '+parameter_file+' --numThreads '+str(args.threads) + ' '
  cmd1 += '--fofn '+tdir+'/rawname.txt -o '+tdir+'/'+output_dir
  sys.stderr.write(cmd1+"\n")
  #os.system(cmd1)
  subprocess.call(cmd1,shell=True)
  return

def setup_temporary_directory(args):
  if args.specific_tempdir:
    tdir = args.specific_tempdir.rstrip('/')
    if not os.path.isdir(args.specific_tempdir):
      os.makedirs(tdir)
    return tdir
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
# post:      An array of dictionaries,
#            with 'name' and 'seq' entires.
# modifies:  none

def read_fasta_into_array(fasta_filename):
  seqs = []
  with open(fasta_filename) as f:
    seq = ''
    head = ''
    prog = re.compile('^>(.+)')
    for line in f:
      line = line.rstrip()
      m = prog.match(line)
      if m:
        if seq != '':
          val = {}
          val['name'] = head
          val['seq'] = seq
          seqs.append(val)
          seq = ''
        head = m.group(1)
      else:
         seq = seq + line
    if seq != '':
      val = {}
      val['name'] = head
      val['seq'] = seq
      seqs.append(val)
  return seqs

if __name__=="__main__":
  main()
