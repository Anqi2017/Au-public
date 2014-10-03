#!/usr/bin/perl -w 
use strict;

# take a fasta piped to standard input
# write a tsv to standard output

my $seq = '';
my $name = '';
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^>([^\t]+)/) {
    my $oldname = $name;
    $name = $1;
    if($seq ne '') {
      print "$oldname\t$seq\n";
    }
    $seq = '';
  } else {
    $line=~s/\s//g;
    $seq .= $line;
  }
}
if($seq ne '') {
  print "$name\t$seq\n";
}
