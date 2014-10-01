#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);
##############
# Reduce a fastq file to a given number of reads
# Input:  fastq file name, output file name, number of reads in output
# Output:  reduced size fastq
# Modifies: File IO

if(scalar(@ARGV) != 3) { die "<input filename> <output filename> <output size>\n"; }
my $infile = shift @ARGV;
my $outfile = shift @ARGV;
my $size = shift @ARGV;
open(INF,$infile) or die;
my $i = 0;
my @inds;
while(my $line1 = <INF>) {
  chomp($line1);
  if($line1=~/^@/) { } else { die "bad header1 format\n"; }
  chomp(my $line2 = <INF>);
  chomp(my $line3 = <INF>);
  if($line3=~/^\+/) { } else { die "bad header2 format\n"; }
  chomp(my $line4 = <INF>);
  push @inds, $i;
  $i++;
}
close INF;
my @rinds = shuffle(@inds);
@inds = ();
my %use;
for(my $j = 0; $j < $size; $j++) {
  $use{$rinds[$j]} = 1;
}
my $k = 0;
open(INF,"$infile") or die;
open(OF,">$outfile") or die;
while(my $line1 = <INF>) {
  chomp($line1);
  if($line1=~/^@/) { } else { die "bad header1 format\n"; }
  chomp(my $line2 = <INF>);
  chomp(my $line3 = <INF>);
  if($line3=~/^\+/) { } else { die "bad header2 format\n"; }
  chomp(my $line4 = <INF>);
  if(exists($use{$k})) {
    print OF "$line1\n$line2\n$line3\n$line4\n";
  }
  $k++;
}
close OF;
