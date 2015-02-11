#!/usr/bin/perl
use strict;

# Read a list of names
# only output fastq entries with those names
if(scalar(@ARGV) != 1) {
  die "cat my.fastq | ./filter_fastq_by_name_list.pl <name list file>\n";
}
my $namefile = shift @ARGV;
open(INF,$namefile) or die;
my %names;
while(my $line = <INF>) {
  chomp($line);
  $names{$line} = 1;
}
close INF;
while(my $line1 = <STDIN>) {
  chomp($line1);
  chomp(my $line2 = <STDIN>);
  chomp(my $line3 = <STDIN>);
  chomp(my $line4 = <STDIN>);
  if($line1=~/^@(.*)$/) {
    my $name = $1;
    if(exists($names{$name})) {
      print "$line1\n$line2\n$line3\n$line4\n";
    }
  }
}

