#!/usr/bin/perl -w
use strict;
if(scalar(@ARGV) !=1) { die "./dp_filter.pl <minimum number of reads on forward and reverse of each allele>\n"; }
my $min = shift @ARGV;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^#/) { print "$line\n"; next; }
  if($line=~/DP4=(\d+),(\d+),(\d+),(\d+)/) {
    if($1 >= $min && $2 >= $min && $3 >= $min && $4 >= $min) {
      print "$line\n";
    }
  }
}
