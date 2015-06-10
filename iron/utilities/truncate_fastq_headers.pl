#!/usr/bin/perl -w
use strict;
# Take as input a fastq file piped to STDIN
# In the output, strip the header after first whitespace
while(my $line1 = <STDIN>) {
  chomp($line1);
  chomp(my $line2 = <STDIN>);
  chomp(my $line3 = <STDIN>);
  chomp(my $line4 = <STDIN>);
  if($line1=~/^(@\S+)/) { $line1 = $1; }
  print "$line1\n$line2\n$line3\n$line4";
}
