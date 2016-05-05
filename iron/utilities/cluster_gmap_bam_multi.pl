#!/usr/bin/perl -w
use strict;
if(scalar(@ARGV) < 5) { die "cluster_gmap_best.pl <launcher> <directory> <output> <index> <maxpaths>\n"; }
my $launcher = shift @ARGV;
my $input = shift @ARGV;
$input=~s/\/$//;
my $output = shift @ARGV;
$output=~s/\/$//;
my $index = shift @ARGV;
my $paths = shift @ARGV;
#if(-d $output) { die "ERROR output alread exists\n"; }
unless(-d $output) { `mkdir $output`; }
chomp(my @files = `ls $input/*.fq`);
my %nums;
foreach my $file (@files) {
  if($file=~/(\d+)\.fq$/) {
    $nums{$1} = 1;
  }
}
chomp(my @ofiles = `ls $output/`);
foreach my $ofile (@ofiles) {
  if($ofile=~/(\d+)\.bam$/) {
    my $num = $1;
    if(exists($nums{$num})) { delete $nums{$num}; }
  }
}
foreach my $num (keys %nums) {
  my $cmd = "$launcher gmap_fasta_to_psl.py --max_paths $paths --bam --threads 16 --gmap_index $index $input/$num.fq $output/$num.bam";
  print "$cmd\n";
}
