#!/usr/bin/perl -w
use strict;
while(my $line = <STDIN>) {
  chomp($line);
  my ($name,$seq) = split(/\t/,$line);
  print ">$name\n$seq\n";
}
