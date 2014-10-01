#!/usr/bin/perl -w
use strict;

##### psl_to_bam.pl #############
# Convert a psl alignment to a bam alignment
# Input: A psl filename
#        A filename for where to write the sorted bam file
#        A filename for a genome index file (something.fa.fai,
#                   can generated with samtools faidx something.fa)
# Output: A bam file.
# Requirements: samtools binaries AND MISC scripts 
#               must be in the path

if(scalar(@ARGV) != 3) { die "<input psl> <output (*.bam)> <indexed genome path (*.fa.fai)>\n"; }

my $fname = shift @ARGV;
my $fout = shift @ARGV;
if($fout=~/^(.*)\.bam$/) { $fout = $1; }
my $genomeindexpath = shift @ARGV;
my $sortcmd = "psl2sam.pl $fname | samtools view -Sb -t ".$genomeindexpath.' - | samtools sort - '.$fout;
`$sortcmd`;
