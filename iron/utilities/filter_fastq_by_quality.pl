#!/usr/bin/perl -w
use strict;
# Pre: 1.  Takes a fastq as input to STDIN
#      2.  The quality score at which or be low you want to remove sequences
#      3.  The number of bases as or below your quality threshold needing to be
#          present before you remove a sequence from the file.
#      4.  Optionally the word "invert" will output the sequences removed by the filter
# Post:       Writes the filtered fastq to STDOUT.  
# Modifies:  None
if(scalar(@ARGV) < 2) { die "cat input.fastq | filter_fastq_based_on_quality.pl  <low quality value to cut on or below> <low quality cutoff count cut on or above> <(optional) the word 'invert' will return the sequences removed by this filter>\n"; }
my $qvalue = shift @ARGV;
my $qcount = shift @ARGV;
my $invert = 0;
if(scalar(@ARGV) == 1) {
  if(lc($ARGV[0]) eq 'invert') {
    $invert = 1;
  }
}
while(my $a1 = <STDIN>) {
  chomp($a1);
  chomp(my $a2=<STDIN>);
  chomp(my $a3=<STDIN>);
  chomp(my $a4=<STDIN>);
  my @chars1 = split(//,$a4);
  my $lowcount1 = 0;
  foreach my $char (@chars1) {
    if(ord($char) <= $qvalue) {
      $lowcount1++;
    }
  }
  if($invert == 1) {
    if($qcount <= $lowcount1) {
      print  "$a1\n$a2\n$a3\n$a4\n";
    } 
  } elsif($invert == 0) {
    if($qcount <= $lowcount1) {
    } else {
      print "$a1\n$a2\n$a3\n$a4\n";
    }
  }
}
