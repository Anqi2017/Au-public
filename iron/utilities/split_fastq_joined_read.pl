#!/usr/bin/perl -w
use strict;

#### split_fastq_joined_read.pl ###########
# It looks like some fastq files dumped from sra
# files have the two sides of paired reads concatonated in one
# fastq file.  if they are the same length and you know what the
# individual lengths should be and the file format is agreeable,
# this will split them into two files.
# Input:  A fastq file name, and a length of one of the two paired ends
# Output: Two fast files same as the basename of the original with
#         _1.fastq or _2.fastq following the base
# Modifies:  File IO

if(scalar(@ARGV) != 2) { die "<filename.fastq> <length of read>\n"; }
my $filename = shift @ARGV;
my $basename = '';
if($filename=~/^(\S+)\.fastq$/) {
  $basename = $1;
} else { die; }
my $length = shift @ARGV;
open(INF,"$filename") or die "could not open $filename\n";

open(OF1,">".$basename."_1.fastq") or die;
open(OF2,">".$basename."_2.fastq") or die;
while(my $line = <INF>) {
  chomp($line);
  my $nline;
  if($line=~/^(@.*\s+length\s*=\s*)(\d+)/) {
    $nline = $1 . $length;
  } else { die "unexpected fastq name line format: $line\n"; }
  print OF1 "$nline\n";
  print OF2 "$nline\n";
  chomp($line = <INF>);
  my $leftseq;
  my $rightseq;
  if($line=~/^(\S{$length,$length})(\S{$length,$length})$/) {
    $leftseq = $1;
    $rightseq = $2;
  } else { die "unexpected fastq sequence line: $line\n"; }
  print OF1 "$leftseq\n";
  print OF2 "$rightseq\n";
  chomp($line=<INF>);
  my $nline2;
  if($line=~/^(\+.*\s+length\s*=\s*)(\d+)/) {
    $nline2 = $1 . $length;
  } else { die "unexpected fastq name2 line format: $line\n"; }
  print OF1 "$nline2\n";
  print OF2 "$nline2\n";
  chomp($line = <INF>);
  my $leftqual;
  my $rightqual;
  if($line=~/^(\S{$length,$length})(\S{$length,$length})$/) {
    $leftqual = $1;
    $rightqual = $2;
  } else { die "unexpected fastq quality line format: $line\n"; }
  print OF1 "$leftqual\n";
  print OF2 "$rightqual\n";
}
close OF1;
close OF2;
