#!/bin/sh
#
#$ -N genetab
#$ -q KA 
#$ -pe smp 1
#$ -cwd

NAME=""

cat ${NAME}1.txt ${NAME}2.txt ${NAME}3.txt ${NAME}4.txt ${NAME}5.txt ${NAME}6.txt ${NAME}7.txt ${NAME}8.txt ${NAME}9.txt ${NAME}10.txt ${NAME}11.txt ${NAME}12.txt ${NAME}13.txt ${NAME}14.txt ${NAME}15.txt ${NAME}16.txt ${NAME}17.txt > ${NAME}all.txt

sort ${NAME}all.txt

uniq ${NAME}all.txt
