#!/usr/bin/perl -w 
use strict;

# take a fastq file or gzipped fastq file
# write a tsv of read name and length of read

if(scalar(@ARGV) != 1) { die "./pacbio_get_longest_fastq_per_molecule.pl myfastq.fq\n"; }
my $filename = shift @ARGV;
my %long;
my %choice;
if($filename=~/\.gz$/) {
  open(INF,"zcat $filename |") or die;
} else {
  open(INF,$filename) or die;
}
while(my $line1 = <INF>) {
  chomp($line1);
  chomp(my $line2=<INF>);
  chomp(my $line3=<INF>);
  chomp(my $line4=<INF>);
  if($line1=~/^@([^\t]+)/) {
    my $name = $1;
    if($line2 ne '') {
      my $pname;
      if($name=~/(^[^\/]+\/\d+)/) {
        $pname = $1;
      } else { die "odd name $name\n"; }
      if(!exists($long{$pname})) {
        $long{$pname} = 0;
        $choice{$pname} = $name;
      }
      if($long{$pname} < length($line2)) {
        $long{$pname} = length($line2);
        $choice{$pname} = $name;
      }
    }
  } 
}
close INF;
if($filename=~/\.gz$/) {
  open(INF,"zcat $filename |") or die;
} else {
  open(INF,$filename) or die;
}
while(my $line1 = <INF>) {
  chomp($line1);
  chomp(my $line2=<INF>);
  chomp(my $line3=<INF>);
  chomp(my $line4=<INF>);
  if($line1=~/^@([^\t]+)/) {
    my $name = $1;
    if($line2 ne '') {
      my $pname;
      if($name=~/(^[^\/]+\/\d+)/) {
        $pname = $1;
      } else { die "odd name $name\n"; }
      if($choice{$pname} eq $name) {
        print "$line1\n$line2\n$line3\n$line4\n";
      }
    }
  } 
}
close INF;
