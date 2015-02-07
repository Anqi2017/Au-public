#!/usr/bin/perl -w 
use strict;
if(scalar(@ARGV) != 1) { die "cat my.sam | ./filter_sam_by_name_list.pl <name list>\n"; }
my $fname = shift @ARGV;
open(INF,$fname);
my %names;
while(my $line = <INF>) {
  chomp($line);
  $names{$line} = 1;
}
close INF;
while(my $line = <STDIN>) {
  chomp($line);
  my @f = split(/\t/,$line);
  if(scalar(@f)<6) {
    print "$line\n";
    next;
  }
  if(exists($names{$f[0]})) {
    print "$line\n";
  }
}
