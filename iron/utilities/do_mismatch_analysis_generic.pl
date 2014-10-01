#!/usr/bin/perl -w
use strict;
use ErrorFlankAnalysis;
use ErrorCodec;
use List::Util qw(shuffle);
if(scalar(@ARGV) != 2) { die; }
my $infile = shift @ARGV;
my $outfile = shift @ARGV;
my @error_codes;

open(INF,$infile) or die;
while(my $line = <INF>) {
  chomp($line);
  my @fields = split(/\t/,$line);
  my $ec = new ErrorCodec($fields[3],$fields[2]); # ccs then sub read
  my $code = $ec->get_error_code();
  push @error_codes, $code;
}
my $efa = new ErrorFlankAnalysis(\@error_codes);
$efa->print_overall_error();
open(OF,">$outfile") or die;
my $r = $efa->get_mismatch_primary_error_report();
print OF $r;
my $r1 = $efa->get_mismatch_flanking_error_report(1);
print OF $r1;
close OF;
