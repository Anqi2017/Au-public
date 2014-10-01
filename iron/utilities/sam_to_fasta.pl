#!/usr/bin/perl -w
use strict;
my $i = 0;
while(my $line = <STDIN>) {
  chomp($line);
  my @f = split(/\t/,$line);
  if($f[9]) {
  if($f[9]=~/[ACTGactguU]/) {
    $i++;
    print ">$i\n$f[9]\n";
  }
  }
}
