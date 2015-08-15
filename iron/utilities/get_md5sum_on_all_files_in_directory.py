#!/usr/bin/python
import sys, argparse, os, subprocess, re, multiprocessing

def main():
  parser = argparse.ArgumentParser(description="for a directory with a ton of files, get checksums")
  parser.add_argument('directory',help="DIRECTORY to check files of")
  parser.add_argument('--threads',type=int,help="INT default cpu_count")
  parser.add_argument('-o','--output',help="FILENAME output file otherwise STDOUT")
  args = parser.parse_args()
  cpus = multiprocessing.cpu_count()
  if args.threads:
    cpus = args.threads
  if not os.path.isdir(args.directory):
    sys.stderr.write("ERROR not a dictory\n")
    return
  files = [f for f in os.listdir(args.directory) if os.path.isfile('/'.join([args.directory.rstrip('/'),f]))]
  pool = multiprocessing.Pool(cpus)
  results = {}
  i = 0
  tot = len(files)
  for file in files:
    loc = '/'.join([args.directory.rstrip('/'),file])
    results[i] = pool.apply_async(calculate_checksum,[loc,i+1,tot])
    i+=1
  pool.close()
  pool.join()
  sys.stderr.write("\n")
  #nums = sorted(results.keys())
  onum = results.keys()
  of = sys.stdout
  if args.output: of = open(args.output,'w')
  for i in sorted(onum):
    of.write(results[i].get()+"\n")
  of.close()
    
def calculate_checksum(loc,z,tot):
  cmd = "md5sum "+loc
  ps = subprocess.Popen(cmd.split(),stdout=subprocess.PIPE)
  output = ps.communicate()[0].rstrip()
  m = re.match('^(\S+)',output)
  sys.stderr.write("\r"+str(z)+'/'+str(tot)+" ")
  return m.group(1) + "\t" + loc

if __name__=="__main__":
  main()
