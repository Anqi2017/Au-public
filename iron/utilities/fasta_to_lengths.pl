#!/usr/bin/perl -w 
use strict;

# take a fasta piped to standard input
# write a tsv of read name and length of read

my $seq = '';
my $name = '';
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^>([^\t]+)/) {
    my $oldname = $name;
    $name = $1;
    if($seq ne '') {
      print "$oldname\t".length($seq)."\n";
    }
    $seq = '';
  } else {
    $line=~s/\s//g;
    $seq .= $line;
  }
}
if($seq ne '') {
  print "$name\t".length($seq)."\n";
}
