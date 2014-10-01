#!/usr/bin/perl -w
use strict;
use SequenceBasics qw(rc);
if(scalar @ARGV < 1) { 
  die "reverse_compliment.pl <seq>\n";
}
my $seq = shift @ARGV;
print rc($seq) . "\n";
