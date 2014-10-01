#!/usr/bin/python

import sys
import os

col_str = "0,0,0"
index = 2
if len(sys.argv) >= 2:
    genephed_filename = sys.argv[1]
    if len(sys.argv)>=3:
        if sys.argv[2].upper() == "B":
            col_str = "0,0,255"
            index = 0
        elif sys.argv[2].upper() == "G":
            col_str = "0,255,0"
            index = 2
        elif sys.argv[2].upper() == "R":
            col_str = "255,0,0"
            index = 1
        elif sys.argv[2].upper() == "Y":
            col_str = "255,255,0"
            index = 2
        elif sys.argv[2].upper() == "P":
            col_str = "190,50,200"
            index = 0

else:
    print("usage:genephed2psl.py genephed_file ")
    print("or gpd2psl.py ")
    sys.exit(1)
################################################################################

#3766    0       0       0       0       0       3       21276   -       chr19:53611131-53636173.1       3766    0       3766    chr19   59128983        53611131        53636173        4       2030,1224,340,172,      0,2030,3254,3594,    53611131,53618462,53625656,53636001,

#chr19:53611131-53636173 chr19:53611131-53636173.1       chr19   -       53611131        53636173        53611131        53636173        4       53611131,53618462,53625656,53636001,    53613161,53619686,53625996,53636173,

#track   name=junctions  description="SpliceMap junctions" itemRgb="On"
#chr19   53611131        53636173        chr19:53611131-53636173.1       50      -       53611131        53636173        255,192,192     4       2030,1224,340,172,      0,7331,14525,24870,

################################################################################
def design_col(score,range_start,range_end,col_str,index):
    ls = col_str.split(',')
    ls[index] = str(int( range_end - (range_end - range_start) * float(score)/1000 ))
    return ','.join(ls)
################################################################################

genephed = open(genephed_filename,'r')

print "track	name=junctions	description=\"SpliceMap junctions\"	itemRgb=\"On\""
for line in genephed:
    ls = line.strip().split("\t")
    N_exon = ls[8]
    exon_start_ls = ls[9].strip(',').split(',')
    exon_end_ls = ls[10].strip(',').split(',')
    block_size_ls = []
    T_start_ls = []
    i = 0
    for start in exon_start_ls:
        end = exon_end_ls[i]
        block_size = int(end)-int(start)
        block_size_ls.append(str(block_size))
        T_start_ls.append(  str( int(start)- int(exon_start_ls[0]) )  )
        i += 1
    
    output_ls = [ls[2]]
    output_ls.append( ls[4] )
    output_ls.append( ls[5] )
    output_ls.append( ls[1] )
    score = "1000"
    output_ls.append( score )
    output_ls.append( ls[3] )
    output_ls.append( ls[4] )
    output_ls.append( ls[5] )
    color = design_col(score, 0,255,col_str,index)
    output_ls.append( color )
    output_ls.append( ls[8] )

    output_ls.append( ','.join(block_size_ls) + ',' )

    output_ls.append( ','.join(T_start_ls) + ',' )

    print '\t'.join(output_ls)
genephed.close()

################################################################################

