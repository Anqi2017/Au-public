#!/usr/bin/perl
use strict;

# Read a list of names
# only output fastq entries with those names
# if you add the optional 'inv' to then it will print out everything except
#    what is on the list.
if(scalar(@ARGV) != 1 && scalar(@ARGV) != 2) {
  die "cat my.fastq | ./filter_fastq_by_name_list.pl <name list file> inv(optional)\n";
}
my $namefile = shift @ARGV;
my $inv = '';
if(scalar(@ARGV) > 0) { $inv = shift @ARGV; }
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
    if($inv ne 'inv') {
      if(exists($names{$name})) {
        print "$line1\n$line2\n$line3\n$line4\n";
      }
    } else {
      if(!exists($names{$name})) {
        print "$line1\n$line2\n$line3\n$line4\n";
      }
    }
  }
}

