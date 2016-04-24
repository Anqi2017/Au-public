#!/usr/bin/perl -w 
use strict;
my $cnt = 0;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^>/) {
    $cnt++;
    $line = ">s_".$cnt;
  }
  print "$line\n";
}
