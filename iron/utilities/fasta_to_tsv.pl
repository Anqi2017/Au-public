#!/usr/bin/perl -w 
use strict;
my $seq = '';
my $name = '';
if(scalar(@ARGV) != 1) { die; }
my $filename = shift @ARGV;
open(OF,">$filename.table") or die;
open(INF,"$filename") or die;
while(my $line = <INF>) {
  chomp($line);
  if($line=~/^>([^\t]+)/) {
    my $oldname = $name;
    $name = $1;
    if($seq ne '') {
      print OF "$oldname\t$seq\n";
    }
    $seq = '';
  } else {
    $line=~s/\s//g;
    $seq .= uc($line);
  }
}
if($seq ne '') {
  print OF "$name\t$seq\n";
}
close INF;
close OF;
