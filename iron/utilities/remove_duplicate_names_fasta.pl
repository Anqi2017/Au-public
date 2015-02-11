#!/usr/bin/perl -w
use strict;
use SequenceBasics qw(read_fasta_into_hash);
my $file = shift @ARGV;
my $f = read_fasta_into_hash($file);
foreach my $name (keys %{$f}) {
  print ">$name\n";
  print $f->{$name}."\n";
}
