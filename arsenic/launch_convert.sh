#Original by Benjamin Deonovic
#Modified by Anhony Rhoads

#!/bin/sh
#
#$ -q UI
#$ -r y
#$ -pe smp 1
#$ -cwd

cat $file | bamToBed -i - > ${file%.bam}.bed

/Shared/Au/anthony/source/CisGenome/cisgenome_project/bin/file_bed2aln -i $file -o $file.aln
