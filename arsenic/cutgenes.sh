#!/bin/sh
#
#$ -N genecut
#$ -q KA 
#$ -pe smp 1
#$ -cwd

NAME=""
NAME2=""

cat ${NAME} | cut -f 1 > /Shared/Au/anthony/wkdir/houston/refseq/tables/${NAME2}


