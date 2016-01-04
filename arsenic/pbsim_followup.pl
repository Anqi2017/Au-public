#!/usr/bin/perl -w
use strict;
##############

#Assign unique names to fastq entries and obtain number of reads per transcript after PBSIM

FASTQS=$(ls | grep fastq)

for f in $FASTQS
do
python /Shared/Au/anthony/utilities/change_readnames.py $f 1 > $f.uniqname
done


MAFS=$(ls | grep maf)
for m in $MAFS
do
cat $m | grep ^a | wc -l > $m.hits
done
