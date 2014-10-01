#!/usr/bin/perl -w
use strict;
############
# What if a fastq read is not long enough?
# Add some N's?
# Why not.
# Pre:  fastq file piped to STDIN with four line pattern
# @readname
# sequence
# +
# pattern
# length of read as argument
# Post: lengthen the read.
if(scalar(@ARGV) < 1) { die; }
my $len = shift @ARGV;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^(@\S+)/) { 
    my $name = $1;
    print "$name\n";
    chomp(my $nline=<STDIN>);
    my $add = 0;
    $add = $len-length($nline);
    if($add < 0) { die "You have a read of length " . length($nline). " which is longer than your desired length of $len.  The desired lenght needs to be longer.\n"; }
    for(my $i = 0; $i < $add; $i++) {
      $nline .= 'N';
    }
    print "$nline\n";
    chomp(my $nline2=<STDIN>);
    if($nline2=~/^(\+)/) { 
      print "$1\n";
    } else { die "format problem, expecting + at beginning of line sorry\n"; }
    chomp(my $nline3=<STDIN>);
    for(my $i = 0; $i < $add; $i++) {
      $nline3 .= 'B';
    }
    print "$nline3\n";
  } 
}
