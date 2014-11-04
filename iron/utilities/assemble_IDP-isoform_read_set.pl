#!/usr/bin/perl -w
use strict;
use SequenceBasics qw(read_fasta_into_hash);
# Pre:
#      CCS95 fasta
#      Pre-LSC fasta sub75_ccs90-95
#      LSC-full fasta
#      LSC-swapped fasta
#      LSC-swapped list
if(scalar(@ARGV) != 5) { die "./assemble_fusion_set.pl <ccs95 fasta> <pre-lsc fasta> <lsc full fasta> <lsc swapped fasta> <lsc swapped list>\n"; }
my $ccs95file = shift @ARGV;
my $prelscfile = shift @ARGV;
my $fullfile = shift @ARGV;
my $swappedfastafile = shift @ARGV;
my $swappedlistfile = shift @ARGV;

open(INF,"$swappedlistfile") or die;
my %swaplist;
while(my $line = <INF>) {
  chomp($line);
  my @f = split(/\t/,$line);
  my $base;
  if($f[0]=~/^([^\/]+\/\d+)/) {
    $base = $1;
  } else { die; }
  $swaplist{$base} = 1;
}
close INF;

my $full = read_fasta_into_hash($fullfile);
foreach my $name (keys %{$full}) {
  my $base;
  if($name=~/^([^\/]+\/\d+)/) {
    $base = $1;
  } else { die; }
  if(!exists($swaplist{$base})) {
    print ">$name\n$full->{$name}\n";
  }
}

my $swapped = read_fasta_into_hash($swappedfastafile);
my %bases;
foreach my $name (keys %{$swapped}) {
  my $base;
  if($name=~/^([^\/]+\/\d+)/) {
    $base = $1;
  } else { die; }
  $bases{$base} = 1;
  print ">$name\n$swapped->{$name}\n";
}

my $prelsc = read_fasta_into_hash($prelscfile);
foreach my $name (keys %{$prelsc}) {
  my $base;
  if($name=~/^([^\/]+\/\d+)/) {
    $base = $1;
  } else { die; }
  if(!exists($bases{$base})) {
    print ">$name\n$prelsc->{$name}\n";
  }  
}

my $ccs95 = read_fasta_into_hash($ccs95file);
foreach my $name (keys %{$ccs95}) {
  print ">$name\n$ccs95->{$name}\n";
}
