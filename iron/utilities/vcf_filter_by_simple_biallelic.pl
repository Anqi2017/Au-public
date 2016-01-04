#!/usr/bin/perl -w
use strict;
if(scalar(@ARGV) !=0) { die "cat my.vcf | ./dp_filter.pl \n"; }
my $min = shift @ARGV;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^#/) { print "$line\n"; next; }
  my @fields = split("\t",$line);
  if($fields[3]=~/^[ACTUG]$/ && $fields[4]=~/^[ACTUG]$/) {
    print "$line\n";
  }
}
