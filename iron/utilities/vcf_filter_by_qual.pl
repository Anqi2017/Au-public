#!/usr/bin/perl -w
use strict;
if(scalar(@ARGV) !=1) { die "./dp_filter.pl <minimum number of reads on forward and reverse of each allele>\n"; }
my $min = shift @ARGV;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^#/) { print "$line\n"; next; }
  my @fields = split(/\t/,$line);
  if($fields[5] >= $min) { print "$line\n";} 
}
