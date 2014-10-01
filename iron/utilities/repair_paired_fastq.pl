#!/usr/bin/perl -w
use strict;
# Pre:  Input files that are mispaired, but share the same name in the first nonwhitespace
# Post: Properly paired files
if(scalar(@ARGV) !=4) { die "<infile1> <infile2> <outfile1> <outfile2>\n"; }
my $infile1 = shift @ARGV;
my $infile2 = shift @ARGV;
my $outfile1 = shift @ARGV;
my $outfile2 = shift @ARGV;
if($infile1=~/\.gz$/) {
  open(INF,"zcat $infile1|") or die;
} else {
  open(INF,$infile1) or die;
}
my %names;
while(my $line1 = <INF>) {
  chomp($line1);
  chomp(my $line2 = <INF>);
  chomp(my $line3 = <INF>);
  chomp(my $line4 = <INF>);
  if(!($line1=~/^@/)) {  die "strange header1 $line1\n"; }
  if(!($line3=~/^\+/)) {  die "strange header2 $line3\n"; }
  my $name;
  if($line1=~/^@(\S+)/) {
    $name = $1;
  }
  my $four = "$line1\n$line2\n$line3\n$line4\n";
  if(exists($names{$name})) { die "problem already saw $name\n"; }
  $names{$name} = $four;
}
close INF;
if($infile2=~/\.gz$/) {
  open(INF,"zcat $infile2|") or die;
} else {
  open(INF,"$infile2") or die;
}
if($outfile1=~/\.gz$/) {
  open(OF1,"| gzip >$outfile1") or die;
} else {
  open(OF1,"$outfile1") or die;
}
if($outfile2=~/\.gz$/) {
  open(OF2,"| gzip >$outfile2") or die;
} else {
  open(OF2,"$outfile2") or die;
}
while(my $line1 = <INF>) {
  chomp($line1);
  chomp(my $line2 = <INF>);
  chomp(my $line3 = <INF>);
  chomp(my $line4 = <INF>);
  if(!($line1=~/^@/)) {  die "strange header1 $line1\n"; }
  if(!($line3=~/^\+/)) {  die "strange header2 $line3\n"; }
  my $name;
  if($line1=~/^@(\S+)/) {
    $name = $1;
  }
  my $four = "$line1\n$line2\n$line3\n$line4\n";
  if(exists($names{$name})) {  
    print OF1 $names{$name};
    print OF2 $four;
  }
}
