#!/usr/bin/perl -w
use strict;
while(my $line = <STDIN>) {
  chomp($line);
  my ($name,$f2,$f3,$f4) = split(/\t/,$line);
  print "@"."$name\n$f2\t$f3\t$f4\n";
}
