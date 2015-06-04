#!/usr/bin/perl -w
use strict;
if(scalar(@ARGV) != 1) { die; }
my $fname = shift @ARGV;
my %chrs;
open(INF,"$fname") or die;
while(my $line = <INF>) {
  chomp($line);
  $chrs{$line} = 1;
}
close INF;
while(my $line = <STDIN>) {
  chomp($line);
  my @f = split(/\t/,$line);
  if(exists($chrs{$f[2]})) {
    print "$line\n";
  }
}  
