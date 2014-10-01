#!/usr/bin/perl -w
use strict;
##############
# Left trim a fastq file by a given number of bases
# Input:  fastq file to STDIN, the word left or right, number of bases to trim
# Output:  trimmed fastq to STDOUT
# Modifies: File IO

if(scalar(@ARGV) != 2) { die "./trim_fastq.pl <left or right> <num to trim>\n"; }
my $leftorright = shift @ARGV;
my $size = shift @ARGV;
if($leftorright ne 'left' && $leftorright ne 'right') { die "left or right please.\n"; }
my $i = 0;
my @inds;
while(my $line1 = <STDIN>) {
  chomp($line1);
  my $len;
  my $head1;
  if($line1=~/^(@.*length\s*=\s*)(\d+)$/) { 
    $head1 = $1;
    $len = $2;
    my $ns = $len - $size;
    $line1=$head1.$ns;
  } elsif ($line1=~/^@/) { # ok if its just a plain header
  } else { die "bad header1 format\n"; }
  chomp(my $line2 = <STDIN>);
  $len = length($line2);
  my $newsize = $len-$size;
  if($len - $size <= 0) { die "error: too short\n"; }
  if($leftorright eq 'left') {
    $line2=substr($line2,-1*$newsize);
  } else {
    $line2=substr($line2,0,$newsize);
  }
  chomp(my $line3 = <STDIN>);
  my $head2;
  if($line3=~/^(\+.*length\s*=\s*)(\d+)$/) { 
    $head2 = $1;
    if($2 != $len) { die "weird fastq.  the length in their header didn't match the length of the sequence\n" };
    $line3 = $head2.$newsize;
  } elsif($line3=~/^\+/) {
  } else { die "bad header2 format\n"; }
  chomp(my $line4 = <STDIN>);
  if($leftorright eq 'left') {
    $line4 = substr($line4,-1*$newsize);
  } else {
    $line4 = substr($line4,0,$newsize);
  }
  print "$line1\n$line2\n$line3\n$line4\n";
}
