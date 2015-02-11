#!/usr/bin/perl -w
use strict;
if(scalar(@ARGV) != 2) { die "./compress_sam.pl <input file> <output file>\n"; }
# Pre:  Arguments input file name (sam), outputfile name (bam)
#       Requires samtools installed in path
# Post: Writes outputfile name bam, a sorted bam file
# Modifiles: uses temp files
my $infile = shift @ARGV;
my $outfile = shift @ARGV;
unless(-d "/tmp/weirathe") { `mkdir /tmp/weirathe`; }
my $rnum = int(rand()*100000000);
my $tfile1 = "/tmp/weirathe/bamA".$rnum.".bam";
my $tfile2base = "/tmp/weirathe/bamB".$rnum;
my $tfile2 = $tfile2base.".bam";
my $cmd1 = "samtools view -bS $infile >  $tfile1";
`$cmd1`;
my $cmd2 = "samtools sort $tfile1 $tfile2base";
`$cmd2`;
my $cmd3 = "mv $tfile2 $outfile";
`$cmd3`;
`rm $tfile1`;
