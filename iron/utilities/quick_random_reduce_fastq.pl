#!/usr/bin/perl -w
use strict;
##############
# Reduce a fastq file to an approximate fraction its current size
# Input:  fraction like 0.0001
# Output:  reduced size fastq
# Modifies: File IO

if(scalar(@ARGV) != 1) { die "<output fraction>\n"; }
my $size = shift @ARGV;
my $i = 0;
my @inds;
while(my $line1 = <STDIN>) {
  chomp($line1);
  if($line1=~/^@/) { } else { die "bad header1 format\n"; }
  chomp(my $line2 = <STDIN>);
  chomp(my $line3 = <STDIN>);
  if($line3=~/^\+/) { } else { die "bad header2 format\n"; }
  chomp(my $line4 = <STDIN>);
  if(rand() <= $size) {
    print "$line1\n$line2\n$line3\n$line4\n";
  }
}
