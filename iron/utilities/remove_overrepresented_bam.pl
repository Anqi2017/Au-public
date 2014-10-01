#!/usr/bin/perl -w
use strict;

###################
# remove overrepresented sequences from a bam
# Pre:  a bam file and choose a threshold for percentage that is acceptable, requires samtools
# Post: a bam file without that sequence
# Modifies: None

my $cnt = 0;
my %seqs;
my $tot = 0;
if(scalar(@ARGV) < 3) { die; }
my $filename = shift @ARGV;
my $outfile = shift @ARGV;
our $threshold = shift @ARGV;
open(INF,"samtools view $filename |") or die;
while(my $line = <INF>) {
  $tot++;
  chomp($line);
  my @fields = split(/\t/,$line);
  if(!exists($seqs{$fields[9]})) {
    $seqs{$fields[9]} = 0;
  }
  $seqs{$fields[9]}++;
}
close INF;

my @over;

my $ntot = get_over(\%seqs,\@over,$tot);

while($ntot < $tot) {
  $tot = $ntot;
  $ntot = get_over(\%seqs,\@over,$tot);
}
my %bad;
foreach my $over (@over) {
  $bad{$over} = 1;
}

open(INF,"samtools view -h $filename |") or die;
open(OF,"| samtools view -bS - > $outfile") or die;
while(my $line = <INF>) {
  chomp($line);
  if($line=~/^@/) {
    print OF "$line\n";
  }
  else {
    my @fields = split(/\t/,$line);
    if(exists($bad{$fields[9]})) {}
    else { print OF "$line\n"; }
  }
}
close OF;
close INF;


sub get_over {
  my $seqref = shift @_;
  my $overref = shift @_;
  my $tot = shift @_;
  foreach my $seq (keys %{$seqref}) {
    my $num = $seqref->{$seq};
    if($num/$tot > $threshold/100) {
      push @{$overref}, $seq;
      $tot-=$num;
      delete $seqref->{$seq};
    }
  }
  return $tot;
}
