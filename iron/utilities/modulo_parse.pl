#!/usr/bin/perl -w
use strict;
##### modulo_parse.pl #####
# Given a dividing factor and an offset, only output those lines.
# Useful for dividing a job up for cluster work quickly.
# Input: dividing factor (example 100 would make 100 files)
#        offset (example 1 will get the first for the first job)
#        And input file input into STDIN
# Output: A subset of the input file in STDOUT
if(scalar(@ARGV) != 2) { die; }
my $n = shift @ARGV;
my $k = shift @ARGV;
my $i = 0;
while(my $line = <STDIN>) {
  $i++;
  if($i % $n==$k) {
    chomp($line);
    print "$line\n";
  }
}
