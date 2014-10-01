#!/usr/bin/perl -w
use strict;
############
# What if a fastq file is not long enough?
# Add some N's?
# Why not.
# Pre:  fastq file piped to STDIN with four line pattern
# @readname
# sequence
# +readname
# pattern
# a number of bases to add as argument
# Post: lengthen the read.
if(scalar(@ARGV) < 1) { die; }
my $add = shift @ARGV;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^(@.*length=\s*)(\d+)(.*)$/) { 
    my ($first, $num, $last) = ($1,$2,$3);
    $num += $add;
    $line = $first.$num.$last;
    print "$line\n"; 
    chomp(my $nline=<STDIN>);
    for(my $i = 0; $i < $add; $i++) {
      $nline .= 'N';
    }
    print "$nline\n";
  } 
  elsif($line=~/^(\+.*length=\s*)(\d+)(.*)$/) { 
    my ($first, $num, $last) = ($1,$2,$3);
    $num += $add;
    $line = $first.$num.$last;
    print "$line\n"; 
    chomp(my $nline=<STDIN>);
    for(my $i = 0; $i < $add; $i++) {
      $nline .= 'B';
    }
    print "$nline\n";
  } 
}
