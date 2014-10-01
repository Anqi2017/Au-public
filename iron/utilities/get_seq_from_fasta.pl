#!/usr/bin/perl -w
use strict;
use SequenceBasics qw(read_fasta_into_hash rc);
if(scalar @ARGV < 1) { 
  die "get_seq_from_fasta.pl <fasta> <chrom> <start optional> <stop optional> <direction optional>"; 
}
my $infile = shift @ARGV;
my $f = read_fasta_into_hash($infile);
if(scalar @ARGV < 1) {
  print "Select a chromosome from the following:\n";
  my $o = '';
  foreach my $chrom (sort {$a cmp $b} keys %{$f}) {
    $o.= "$chrom, ";
  }
  chop($o);
  chop($o);
  print "$o\n";
  die "get_seq_from_fasta.pl <fasta> <chrom> <start optional> <stop optional> <direction optional>"; 
}
my $chrom = shift @ARGV;
if(!exists($f->{$chrom})) {
  print "Select a chromosome from the following:\n";
  my $o = '';
  foreach my $chrom (sort {$a cmp $b} keys %{$f}) {
    $o.= "$chrom, ";
  }
  chop($o);
  chop($o);
  print "$o\n";
  die "get_seq_from_fasta.pl <fasta> <chrom> <start optional> <stop optional>"; 
}
if(scalar(@ARGV) == 0) { # print whole chromosome
  print $f->{$chrom} . "\n";
} elsif(scalar(@ARGV) == 1) {
  my $start = shift @ARGV;
  print substr($f->{$chrom},$start-1,1) . "\n";
} elsif(scalar(@ARGV) >= 2) {
  my $start = shift @ARGV;
  my $stop = shift @ARGV;
  my $dir = '+';
  if(scalar(@ARGV) == 1) {
    $dir = shift @ARGV;
  }
  my $seq = substr($f->{$chrom},$start-1,$stop-$start+1);
  if($dir eq '-') { $seq = rc($seq); }
  print $seq . "\n";
} else {
  die "get_seq_from_fasta.pl <fasta> <chrom> <start optional> <stop optional>"; 
}


