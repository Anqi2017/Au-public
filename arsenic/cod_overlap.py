#!/usr/bin/python
import os
import sys
from random import randint
import re
import subprocess
#Input: two .cod files, a cutoff rate 
#       Requires bedtools installed to path
#       Requires write access to /tmp/directory
#Output: intersection of the two .cod files in bed format 
#        Output format is tsv:
#        <file1 chr> <file1 start> <file1 stop> <file2 chr> <file2 start> <file2 stop> <overlap (base paircount)>
#        coordinates are 1-indexed
#Modifies: creates a temporary directory inside of /tmp/ ie /tmp/ajrhoads
#          creates temporary files in that directory, deletes them when finished
#          writes to STDOUT



def main():
    if not os.path.exists("/tmp/ajrhoads"):
        os.makedirs("/tmp/ajrhoads")

    if len(sys.argv) != 4:
        print sys.argv[0]+" <codfile 1> <codfile 2> <overlap fraction (i.e. 0.9)>"
        return 

    fileACod=sys.argv[1]
    fileBCod=sys.argv[2]
    fractionOverlap=sys.argv[3]

    fileABed="/tmp/ajrhoads/overlap1_"+str(randint(1,100000000))+".bed"
    fileBBed="/tmp/ajrhoads/overlap2_"+str(randint(1,100000000))+".bed"

    codtobed(fileACod, fileABed)
    codtobed(fileBCod, fileBBed)

    cmd = 'bedtools intersect -a '+fileABed+' -b '+fileBBed+' -wo -f '+fractionOverlap+' -r'
    p = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    out, err = p.communicate()

    for line in out.split("\n"):
        line = line.rstrip()
        if re.match('^$', line):
            continue
        print line

    os.remove(fileABed)
    os.remove(fileBBed)




def codtobed(fileNameCod, fileNameBed):
    of=open(fileNameBed, 'w')
    
    with open(fileNameCod) as infile:
        for line in infile:
            if re.match('^#',line):
                continue

            f=line.rstrip().split("\t")

            chr=f[1]
            start=int(f[2])+1
            end=int(f[3])+1

#            if end < start:
#                print "negative strand"

            of.write(chr+"\t"+str(start)+"\t"+str(end)+"\n")

    of.close()

main()


