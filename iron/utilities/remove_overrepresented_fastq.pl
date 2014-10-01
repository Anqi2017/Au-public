#!/usr/bin/perl -w
use strict;
#use Digest::SHA1 qw(sha1_hex);

###################
# remove overrepresented sequences from a fastq
# Pre:  a fastq file and choose a threshold for percentage that is acceptable, requires samtools
# Post: a fastq file without that sequence
# Modifies: None

my $cnt = 0;
my %seqs;
my $tot = 0;
if(scalar(@ARGV) < 4) { die "<infile> <outfile> <threshold> <length to look>\n"; }
my $filename = shift @ARGV;
my $outfile = shift @ARGV;
our $threshold = shift @ARGV;
our $look = shift @ARGV;
open(INF,"$filename") or die;
while(my $line1 = <INF>) {
  $tot++;
  chomp($line1);
  chomp(my $line2 = <INF>);
  #my $val = sha1_hex(substr($line2,0,$look));
  my $val = substr($line2,0,$look);
  if(!exists($seqs{$val})) {
    $seqs{$val} = 0;
  }
  $seqs{$val}++;
  chomp(my $line3 = <INF>);
  chomp(my $line4 = <INF>);
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

open(INF,"$filename") or die;
open(OF,">$outfile") or die;
open(OFB,">$outfile.bad") or die;
while(my $line1 = <INF>) {
  chomp($line1);
  chomp(my $line2 = <INF>);
  chomp(my $line3 = <INF>);
  chomp(my $line4 = <INF>);
  my $set = "$line1\n$line2\n$line3\n$line4";
  #my $val = sha1_hex(substr($line2,0,$look));
  my $val = substr($line2,0,$look);
  if(exists($bad{$val})) { print OFB "$set\n"; }
  else { print OF "$set\n"; }
}
close OF;
close OFB;
close INF;


sub get_over {
  my $seqref = shift @_;
  my $overref = shift @_;
  my $tot = shift @_;
  foreach my $seq (keys %{$seqref}) {
    my $num = $seqref->{$seq};
    #print "$seq\n$num\t$tot\n\n";
    if($num/$tot > $threshold/100) {
      push @{$overref}, $seq;
      $tot-=$num;
      delete $seqref->{$seq};
    }
  }
  return $tot;
}
