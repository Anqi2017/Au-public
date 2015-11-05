#!/usr/bin/perl -w
use strict;
#Pre: A psl file streamed in STDIN
#Post: A psl file ordered by the sequence name streamed to STDOUT

open(STREAM,"| sort -k 10,10 > /dev/stdout") or die;
my $z = 0;
while(my $line = <STDIN>) {
  $z+=1;
  if($z==1 && ($line=~/^#/ || !($line=~/^\d+/))) { 
    print $line;
    next; #Skip comment or header
  } 
  print STREAM $line;
}
