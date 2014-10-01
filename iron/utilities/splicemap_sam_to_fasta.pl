#!/usr/bin/perl -w
use strict;
my $i = 0;
while(my $line = <STDIN>) {
  chomp($line);
  my @f = split(/\t/,$line);
  if($f[9]) {
  if($f[9]=~/[ACTGactguU]/) {
    $i++;
    my $name = $i;
    if($f[0]=~/^([^\[]+)/) { $name = $1; }
    print ">$name\n$f[9]\n";
  }
  }
}
