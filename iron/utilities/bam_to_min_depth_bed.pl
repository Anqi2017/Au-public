#!/usr/bin/perl -w 

#########
# Filter a bam file based on the depth
# Input: input bam file name, output bed file name, depth
# Output:  write out the the bed file of locations with depth 10
# Modifies: file io, makes temp files, stderr printing too

use strict;
my $bamfile = shift @ARGV;
my $outbed = shift @ARGV;
my $depth = shift @ARGV;
open(INF,"bedtools genomecov -bg -ibam $bamfile|") or die;
open(OF,"|bedtools merge -i - > $outbed") or die;
while(my $line = <INF>) {
  chomp($line);
  my @fields = split(/\t/,$line);
  if($fields[3] >= 10) {
    print OF "$line\n";
  }
}
close INF;
close OF;
