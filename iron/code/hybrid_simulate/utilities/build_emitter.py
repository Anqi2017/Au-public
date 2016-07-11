#!/usr/bin/python
import argparse, sys, os, gzip, zlib, base64, pickle
from shutil import rmtree
from multiprocessing import cpu_count
from tempfile import mkdtemp, gettempdir
from Bio.Structure import Transcriptome, Transcript
from Bio.Format.Fasta import FastaData
from Bio.Format.GPD import GPD

def main(args):
  sys.stderr.write("Reading reference fasta\n")
  ref_genome = FastaData(open(args.reference_fasta,'rb').read())
  sys.stderr.write("Reading in transcriptome\n")
  output = {}
  txome = Transcriptome()
  z = 0
  with open(args.reference_gpd) as inf:
    for line in inf:
      z+=1
      if z%1000==0:  sys.stderr.write(str(z)+"       \r")
      gpd = GPD(line)
      gpd.set_sequence(ref_genome)
      txome.add_transcript(gpd)
  sys.stderr.write("\n")
  sys.stderr.write("Serializing transcriptome\n")
  output['txome'] = txome.dump_serialized()
  txweights = {}
  weight_type = 'uniform_distribution' #default
  if args.expression_table:
    weight_type = 'expression_table'
    inf = None
    if args.expression_table[-3:]=='.gz':
      inf = gzip.open(args.expression_table)
    else: inf = open(args.expression_table)
    for line in inf:
      f = line.rstrip().split("\t")
      txweights[f[0]] = float(f[1])
  elif args.exponential_distribution: weight_type = 'exponential_distribution'
  output['weight_type'] = weight_type
  output['weights'] = txweights #only matters for expression based
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  of.write(base64.b64encode(zlib.compress(pickle.dumps(output)))+"\n")
  of.close()


  # Temporary working directory step 3 of 3 - Cleanup
  if not args.specific_tempdir:
    rmtree(args.tempdir)

def do_inputs():
  # Setup command line inputs
  parser=argparse.ArgumentParser(description="Establish the transcriptome to emit.  The type of read generated, and perterberations are defined when you emit.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('reference_gpd',help="INPUT GenePred File")
  parser.add_argument('reference_fasta',help="Reference fasta")
  
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--uniform_distribution',default=True,action='store_true',help="Equal probability of any transcript")
  group.add_argument('--exponential_distribution',action='store_true',help="similar to a real distribution")
  group.add_argument('--expression_table',help="TSV table of <transcript name> <expression> to weight expression")

  # Performance details
  parser.add_argument('-o','--output',help="OUTPUTFILE or STDOUT")
  #parser.add_argument('--threads',type=int,default=cpu_count(),help="INT number of threads to run. Default is system cpu count")

  # Temporary working directory step 1 of 3 - Definition
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--tempdir',default=gettempdir(),help="The temporary directory is made and destroyed here.")
  group.add_argument('--specific_tempdir',help="This temporary directory will be used, but will remain after executing.")
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


def external_cmd(cmd,version=None):
  cache_argv = sys.argv
  sys.argv = cmd.split()
  args = do_inputs()
  main(args)
  sys.argv = cache_argv

if __name__=="__main__":
  args = do_inputs()
  main(args)
