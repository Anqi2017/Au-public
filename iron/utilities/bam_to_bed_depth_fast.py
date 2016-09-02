#!/usr/bin/python
import sys, argparse, tempfile, os, StringIO, gzip
from subprocess import PIPE, Popen
from multiprocessing import Pool, cpu_count
from Bio.Format.Sam import BAMFile, SAM

ps = None

def main():
  parser = argparse.ArgumentParser(description="",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input',help="BAM file input")
  parser.add_argument('--threads',type=int,default=cpu_count(),help="Thread count")
  parser.add_argument('--tempdir',help="location of temporary directory to use")
  parser.add_argument('-o','--output',help="Output file name")
  args = parser.parse_args()
  if args.tempdir: args.tempdir = args.tempdir.rstrip('/')
  bf = BAMFile(args.input)
  seqs = bf.get_header().get_sequence_lengths()
  f = tempfile.NamedTemporaryFile(delete=False)
  for seq in seqs:
    f.write(seq+"\t"+str(seqs[seq])+"\n")
  f.close()
  bf.close()
  fout = tempfile.NamedTemporaryFile(delete=False)
  cmd = 'sort -k 1,1 -k2,2n -k3,3n -S4G --parallel='+str(args.threads)
  if args.tempdir: cmd += ' -T '+args.tempdir
  global ps
  ps = Popen(cmd.split(),stdin=PIPE,stdout=fout)

  if args.threads > 1:
    poo = Pool(processes=args.threads)
  for seq in seqs:
    if args.threads > 1:
      poo.apply_async(do_seq,args=(seq,args,f.name),callback=do_output)
    else:
      res = do_seq(seq,args,f.name)
      do_output(res)
  if args.threads > 1:
    poo.close()
    poo.join()
  ps.communicate()
  fout.close()
  of = sys.stdout
  if args.output:
    if args.output[-3:]=='.gz':
      of = gzip.open(args.output,'w')
    else:
      of = open(args.output,'w')
  cmd = 'sort -k 1,1 -k2,2n -k3,3n -S4G --parallel='+str(args.threads)
  inf = open(fout.name)
  p = Popen(cmd.split(),stdout=PIPE,stdin=inf)
  for line in p.stdout:
    of.write(line)
  p.communicate()
  inf.close()
  of.close()
  os.unlink(f.name)
  os.unlink(fout.name)

def do_output(res):
  if not res: return
  global ps
  inf = StringIO.StringIO(res)
  for line in inf:
    ps.stdin.write(line)
  inf.close()

def do_seq(seq,args,fname):
  cmd = 'samtools view '+args.input+' '+seq
  p = Popen(cmd.split(),stdout=PIPE)
  cmd = 'bedtools genomecov -i - -bg -g '+fname
  po2 = Popen(cmd.split(),stdin=PIPE,stdout=PIPE)
  cmd = 'sort -k 1,1 -k2,2n -k3,3n -S1G --parallel='+str(args.threads)
  if args.tempdir: cmd += ' -T '+args.tempdir
  po1 = Popen(cmd.split(),stdin=PIPE,stdout=po2.stdin)
  for line in p.stdout:
    sam = SAM(line)
    if not sam.is_aligned(): continue
    for ex in sam.get_target_transcript(min_intron=68).exons:
      bed = "\t".join([str(x) for x in ex.get_range().get_bed_array()])
      po1.stdin.write(bed+"\n")
  po1.communicate()
  res = po2.communicate()[0]
  return res

if __name__=="__main__":
  main()
