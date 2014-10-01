#!/usr/bin/perl -w
use strict;

use HomopolymerCompression;
if(scalar(@ARGV) != 2) { die "<infile> <outbase>\n"; }
my $infile = shift @ARGV;
my $outbase = shift @ARGV;

my $head = '';
my $seq = '';
open(INF,"$infile") or die;
open(OF,">$outbase") or die;
open(OFLENS,">$outbase.compressionlengths") or die;
open(OFSHORT,">$outbase.compressionshortcuts") or die;
while(my $line = <INF>) {
  chomp($line);
  if($line=~/^>/) {
    if($seq ne '') {
      my $hc = new HomopolymerCompression;
      $hc->load_uncompressed_nucleotides($seq);
      my ($cps,$lens) = $hc->get_compressed_nucleotides();
      print OF "$head\n$cps\n";
      print OFLENS "$head\n$lens\n";
      my $sc = join('|',@{$hc->make_decompression_shortcuts()});
      print OFSHORT "$head\n$sc\n";
    }
    $head = $line;
    $seq = '';
  } else {
    $seq .= $line;
  }
}
if($seq ne '') {
  my $hc = new HomopolymerCompression;
  $hc->load_uncompressed_nucleotides($seq);
  my ($cps,$lens) = $hc->get_compressed_nucleotides();
  print OF "$head\n$cps\n";
  print OFLENS "$head\n$lens\n";
  my $sc = join('|',@{$hc->make_decompression_shortcuts()});
  print OFSHORT "$head\n$sc\n";
}
close OF;
close OFLENS;
close OFSHORT;
close INF;
