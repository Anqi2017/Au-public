#!/usr/bin/python

import sys
import os

if len(sys.argv) >= 2:
    gpd_filename = sys.argv[1]
else:
    print("usage: python remove_alt_chromosomes.py gpd_file > gtf_file_out")
    sys.exit(1)

trans=open(gpd_filename,'r') 

lines = trans.readlines()

for i in range(0,len(lines)):
    line = lines[i]
    if  not (line.startswith('chr1_') or line.startswith('chr2_') or line.startswith('chr3_') or line.startswith('chr4_') or line.startswith('chr5_') or line.startswith('chr6_') or line.startswith('chr7_') or line.startswith('chr8_') or line.startswith('chr9_') or line.startswith('chr10_') or line.startswith('chr11_') or line.startswith('chr12_') or line.startswith('chr13_') or line.startswith('chr14_') or line.startswith('chr15_') or line.startswith('chr16_') or line.startswith('chr17_') or line.startswith('chr18_') or line.startswith('chr19_') or line.startswith('chr20_') or line.startswith('chr21_') or line.startswith('chr22_')  or line.startswith('chrX_') or line.startswith('chrUn') ):
        print(line),

trans.close()
