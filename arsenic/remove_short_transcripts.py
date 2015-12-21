#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 3:
    transcript_filename = sys.argv[1]
    minlength = sys.argv[2]
else:
    print("usage: python remove_short_transcripts.py transcript_file min_length > transcript_file_out")
    sys.exit(1)

trans=open(transcript_filename,'r') 

lines = trans.readlines()

for i in range(0,len(lines)):
    line = lines[i]
    if line.startswith('>'):
        if len(lines[i+1]) > 100: #Not >= because len adds 1
            print(line),
            print(lines[i+1]),

trans.close()
