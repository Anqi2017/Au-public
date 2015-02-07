#!/usr/bin/perl -w 
use strict;

# take a fasta file or gzipped fasta file
# write a tsv of read name and length of read

if(scalar(@ARGV) != 1) { die "./pacbio_get_longest_fasta_per_molecule.pl myfasta.fa\n"; }
my $filename = shift @ARGV;
my $seq = '';
my $name = '';
my %long;
my %choice;
if($filename=~/\.gz$/) {
  open(INF,"zcat $filename |") or die;
} else {
  open(INF,$filename) or die;
}
while(my $line = <INF>) {
  chomp($line);
  if($line=~/^>([^\t]+)/) {
    my $oldname = $name;
    $name = $1;
    if($seq ne '') {
      my $pname;
      if($oldname=~/(^[^\/]+\/\d+)/) {
        $pname = $1;
      } else { die "odd name $oldname\n"; }
      if(!exists($long{$pname})) {
        $long{$pname} = 0;
        $choice{$pname} = $oldname;
      }
      if($long{$pname} < length($seq)) {
        $long{$pname} = length($seq);
        $choice{$pname} = $oldname;
      }
    }
    $seq = '';
  } else {
    $line=~s/\s//g;
    $seq .= $line;
  }
}
if($seq ne '') {
      my $pname;
      if($name=~/(^[^\/]+\/\d+)/) {
        $pname = $1;
      } else { die "odd name $name\n"; }
      if(!exists($long{$pname})) {
        $long{$pname} = 0;
        $choice{$pname} = $name;
      }
      if($long{$pname} < length($seq)) {
        $long{$pname} = length($seq);
        $choice{$pname} = $name;
      }
}
close INF;
$seq = '';
$name = '';
if($filename=~/\.gz$/) {
  open(INF,"zcat $filename |") or die;
} else {
  open(INF,$filename) or die;
}
while(my $line = <INF>) {
  chomp($line);
  if($line=~/^>([^\t]+)/) {
    my $oldname = $name;
    $name = $1;
    if($seq ne '') {
      my $pname;
      if($oldname=~/(^[^\/]+\/\d+)/) {
        $pname = $1;
      } else { die "odd name $oldname\n"; }
      if($choice{$pname} eq $oldname) {
        print ">$oldname\n".$seq."\n";
      }
    }
    $seq = '';
  } else {
    $line=~s/\s//g;
    $seq .= $line;
  }
}
if($seq ne '') {
      my $pname;
      if($name=~/(^[^\/]+\/\d+)/) {
        $pname = $1;
      } else { die "odd name $name\n"; }
      if($choice{$pname} eq $name) {
        print ">$name\n".$seq."\n";
      }
}
close INF;
