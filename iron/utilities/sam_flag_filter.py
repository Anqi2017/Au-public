#!/usr/bin/python

import sys, re, os, inspect

#pre: 1.) a sam file piped to stdin
#     2.) a flag you want to filter on in hex (i.e. 0x100 or 0x800)
#     3.) 1 if you want to print on a flag of 1, or 0 if you want to print on a flag of 0
#post: a sam file

#bring in the folder to the path for our modules
pythonfolder_loc = "../pythonlib"
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0],pythonfolder_loc)))
if cmd_subfolder not in sys.path:
  sys.path.insert(0,cmd_subfolder)
import sam_basics


def main():
  if len(sys.argv) != 3:
    print 'cat myfile.sam | ./'+sys.argv[0] + ' <hex flag (i.e. 0x800)> <integer 1 for print on 1, 0 for print on 0>'
    sys.exit()
  flag = int(sys.argv[1],16)
  todobool = int(sys.argv[2])
  for line in sys.stdin:
    line = line.rstrip()
    if(sam_basics.is_header(line)):
      print line
      continue
    sam = sam_basics.sam_line_to_dictionary(line)
    flagresult = sam_basics.check_flag(sam['flag'],flag)
    if flagresult and todobool == 1:
      print line
    elif not flagresult and todobool == 0:
      print line
main()
