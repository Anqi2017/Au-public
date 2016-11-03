#!/usr/bin/python
import argparse, os, sys, re, multiprocessing
import ONTBasics

z = 0
z2d =0
ztemplate = 0
zcomplement =0
zfail = 0
ztotal = 0
of = None

def main():
  parser = argparse.ArgumentParser(description="Extract fastq from ONT data")
  parser.add_argument('input_file_or_directory',help='FAST5 file with read information or a directory containing fast5 files')
  parser.add_argument('-o','--output',help='FILENAME specify an output FASTQ file')
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--get_2D',action='store_true')
  group.add_argument('--get_tem',action='store_true')
  group.add_argument('--get_com',action='store_true')
  group1 = parser.add_mutually_exclusive_group(required=True)
  #group1.add_argument('--trim_header',action='store_true',help="BOOL True if you want to make sure the name of each read does not contain any spaces (truncates after first whitespace)")
  group1.add_argument('--first_name',action='store_true')
  group1.add_argument('--second_name',action='store_true')
  parser.add_argument('--append_suffix',choices=['pass','fail'],help="you have to set the call but it will add 2D com or tem")
  parser.add_argument('--threads',type=int,default=multiprocessing.cpu_count(),help="defaults to cpu_count")
  args = parser.parse_args()

  global of
  of = sys.stdout
  if args.output: of = open(args.output,'w')

  files = []
  if os.path.isfile(args.input_file_or_directory):
    files.append(args.input_file_or_directory)
    sys.stderr.write("Performing Fastq extraction on one file\n")
  elif os.path.isdir(args.input_file_or_directory):
    sys.stderr.write("Performing Fastq extraction on one directory\n")
    files = ['/'.join([args.input_file_or_directory.rstrip('/'),f]) for f in os.listdir(args.input_file_or_directory) if os.path.isfile('/'.join([args.input_file_or_directory.rstrip('/'),f]))]
  else: 
    sys.stderr.write("input is neither a file nor a directory\n")
    sys.exit()
  fast5_files = []
  for file in files:
    if re.search('\.fast5$',file): fast5_files.append(file)
  if len(fast5_files) == 0:
    sys.stderr.write("no fast5 files given as input\n")
    sys.exit()

  if args.threads > 1:
    cpus = multiprocessing.cpu_count()
    p = multiprocessing.Pool(processes=cpus)
  global ztotal
  ztotal = len(fast5_files)
  for filename in fast5_files:
    if args.threads > 1:
      p.apply_async(do_file,args=(filename,args,),callback=collect_results)
    else:
      r = do_file(filename,args)
      collect_results(r)
    #collect_results(do_file(filename,args))
  if args.threads > 1:
    p.close()
    p.join()
  sys.stderr.write("\n")

def do_file(filename,args):
    y2d = 0
    ytemplate = 0
    ycomplement = 0
    yfail = 0
    obf = ONTBasics.fast5(filename)
    #if args.trim_header:
    #  obf.set_trim_header_bool(True)
    read = None
    extra = ''
    if args.append_suffix:
      extra = '_'+args.append_suffix+'_'
    if obf.extract_2D() and args.get_2D:
      read = obf.extract_2D().fastq()
      extra += '2D'
      y2d=1
    elif obf.extract_template() and args.get_tem:
      read = obf.extract_template().fastq()
      extra += 'tem'
      ytemplate = 1
    elif obf.extract_complement() and args.get_com:
      read = obf.extract_complement().fastq()
      extra += 'com'
      ycomplement = 1
    else:
      #sys.stderr.write("Unable to extract read from "+filename+"\n")
      yfail = 1
    #sys.stderr.write("Extracted file "+str(z)+"/"+str(ztotal)+" ")
    #sys.stderr.write("2D: "+str(z2d)+"  Temp: "+str(ztemplate) + "  Comp: "+str(zcomplement)+"  Fail: "+str(zfail)+"\r")
    obf.close()
    return [read,y2d,ytemplate,ycomplement,yfail,extra,args]

def collect_results(result):
  [read,y2d,ytemplate,ycomplement,yfail,extra,args] = result
  global z
  z += 1
  global z2d
  z2d += y2d
  global ztemplate
  ztemplate += ytemplate
  global zcomplement
  zcomplement += ycomplement
  global zfail
  zfail += yfail
  global ztotal
  sys.stderr.write("Extracted file "+str(z)+"/"+str(ztotal)+" ")
  sys.stderr.write("2D: "+str(z2d)+"  Temp: "+str(ztemplate) + "  Comp: "+str(zcomplement)+"  Fail: "+str(zfail)+"\r")
  global of
  if read:
    lines = read.rstrip("\n").split("\n")
    m = re.match('@(\S+)\s+(\S+)',lines[0])
    name = m.group(1)+extra
    if args.second_name:
      name=m.group(2)+extra
    lines[0]='@'+name
    of.write("\n".join(lines)+"\n")
    #of.write(read)

if __name__=="__main__":
  main()
