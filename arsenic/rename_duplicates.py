#!/usr/bin/python

from __future__ import print_function
import sys
import os

#Adds Look for duplicate transcripts and rename them

if len(sys.argv) >= 3:
    fastq_filename = sys.argv[1]
    duplicates_file = sys.argv[2]
else:
    print("usage: python rename_duplicates.py fastq_filename duplicates_file > fastq_out")
    sys.exit(1)

reads=open(fastq_filename,'r') 
duplicates=open(duplicates_file,'r')

lines = reads.readlines()
dups = duplicates.readlines()


#print duplicates as separate entries

for d in range(0,len(dups)):
  entries = 1
  for i in range(0,len(lines)):
    line = lines[i]
    if i % 2 == 0:

      if line == dups[d]:	
        print (line.rstrip(),'_',entries,sep='')
	entries = entries + 1

	i=i+1
	line=lines[i]
	print (line.rstrip()),


#print non-duplicates

for i in range(0,len(lines)):
  line=lines[i]
  dup=0 # number of duplicates of this line

  for d in range(0,len(dups)):
    if i % 2 == 0:

      if line==dups[d]:
        dup=dup+1

  if i % 2 == 0:
    if dup==0:
      print (line.rstrip()),
      
      i=i+1
      line=lines[i]
      print (line.rstrip()),


reads.close()
duplicates.close()
