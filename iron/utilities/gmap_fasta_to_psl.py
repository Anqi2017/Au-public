#!/usr/bin/python
import argparse, sys, multiprocessing, random, os, re, subprocess
from shutil import rmtree

def main():
  parser = argparse.ArgumentParser(description="Launch GMAP and run it in a temp directory until output is finished.")
  parser.add_argument('input_fasta',help="FASTA file or - for STDIN")
  parser.add_argument('output_psl')
  parser.add_argument('--best',action='store_true',help="Only output the best path")
  parser.add_argument('--gmap_index',required=True,help="Path to gmap index (directory)")
  parser.add_argument('--threads',type=int,default=multiprocessing.cpu_count())
  parser.add_argument('--max_paths',type=int,help="Maximum number of paths to show.")
  parser.add_argument('--max_intron_length',type=int,help="Maximum length of intron.")
  parser.add_argument('--tempdir',default='/tmp')
  args = parser.parse_args()

  args.gmap_index = args.gmap_index.rstrip('/')
  m = re.match('^(.*)$',args.gmap_index)
  if m:
    gmapindexpath = m.group(1)
  m = re.search('([^\/]+)$',args.gmap_index)
  if m:
    gmapindexname = m.group(1)
  else:
    sys.stderr.write("problem reading gmap index "+args.gmap_index+"\n")
    sys.exit()

  maxpathpart = ''
  if args.max_paths:
    maxpathpart = ' -n '+str(args.max_paths)+' '

  maxintronpart = ''
  if args.max_intron_length:
    maxintronpart = ' -K '+str(args.max_intron_length)+' '


  rnum = random.randint(1,10000000)

  args.tempdir = args.tempdir.rstrip('/')+'/weirathe.'+str(rnum)
  if not os.path.exists(args.tempdir):
    os.makedirs(args.tempdir)
  if args.input_fasta == '-':
    args.input_fasta = args.tempdir+'/input.fasta'
    of = open(args.tempdir+'/input.fasta','w')
    for line in sys.stdin:
      of.write(line)
    of.close()
  gmap_cmd = 'gmap --ordered -D '+gmapindexpath+' -f 1 -d '+gmapindexname+' -t '+str(args.threads)+' '+maxpathpart+maxintronpart+' '+args.input_fasta
  sys.stderr.write("executing:\n"+gmap_cmd+"\n")

  allpsl = args.tempdir+'/all.psl'
  #subprocess.call(gmap_cmd+' > '+allpsl,shell=True,stderr=subprocess.PIPE)
  ofe = open('/dev/null','w')
  subprocess.call(gmap_cmd+' > '+allpsl,shell=True,stderr=ofe)

  if args.best:
    vals = {}
    with open(allpsl) as inf:
      for line in inf:
        if re.match('^#',line): 
          continue
        f = line.rstrip().split("\t")
        if f[9] not in vals:
          vals[f[9]]=0
        if int(f[0]) > vals[f[9]]: vals[f[9]] = int(f[0])
          
  of = sys.stdout
  if args.output_psl != '-':
    of = open(args.output_psl,'w')
  with open(allpsl) as inf:
    for line in inf:
      f = line.rstrip().split("\t")
      if args.best:
        if f[9] in vals:
          if int(f[0]) >= vals[f[9]]:
            del vals[f[9]]
            of.write(line)
      else:
        of.write(line)
  of.close()
  rmtree(args.tempdir)

if __name__=='__main__':
  main()
