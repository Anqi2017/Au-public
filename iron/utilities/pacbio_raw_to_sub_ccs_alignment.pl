#!/usr/bin/perl -w
use strict;
use BlasrInterface;
use threads;

###### get_to_ccs_alignments_from_h5.pl ###########
# Use stand alone blasr to find the alignments of subreads to their
# respective ccs reads.  This is a little wonky because as coded we
# align each to the reference genome first and then use the
# reference genome as a guide to align the subread to the ccs for
# the longest ccs hit to the reference genome where it overlaps with
# the subread.
#
# Input: *.bax.h5 filename, reference fasta filename, suffixarray produced by reference fasta file (using sawriter from blasr package), output filename.
# Output: A TSV report that contains the alignment of the sub read to the ccs read, the subread orientation is kept the same as it was in the h5 file.
# Modifies: temporary files in /tmp/

if(scalar(@ARGV) != 4) { 
  my $deathstring =  "<*.bax.h5> <reference fasta> <reference suffix array> <output file>\n"; 
  $deathstring .= "output report is in the following format:\n";
  $deathstring .= "read  errorrate  subread  ccsread  match_count  mismatch_count  numSubreads\n";
  die $deathstring;
}

my $h5 = shift @ARGV;
my $reference = shift @ARGV;
my $sa = shift @ARGV;
my $outreport = shift @ARGV;

my $rnum = int(rand()*10000000);
my $ccsoutput = "/tmp/weirathe.ccs.$rnum.aln";
my $suboutput = "/tmp/weirathe.sub.$rnum.aln";

my $bi = new BlasrInterface($h5,$reference,$sa);

#writes alignments to these files, can do simulataneously
my $thr1 = threads->create(\&write_sub,$bi,$suboutput);
my $thr2 = threads->create(\&write_ccs,$bi,$ccsoutput);
$thr1->join(); #wait for them both to finish before moving on
$thr2->join();


$bi->combine_ccs_and_subreads($ccsoutput,$suboutput);

#get rid of the reference alignment, and align sub to ccs
my @report = $bi->align_subs_to_best_ccs_hits_report();
open (OF,">$outreport") or die;
foreach my $line (@report) {
  chomp($line);
  print OF "$line\n";
}
close OF;

`rm $ccsoutput`;
`rm $suboutput`;


sub write_sub {
  my $bi = shift @_;
  my $suboutput = shift @_;
  $bi->get_blasr_subread_to_genome_alignments_from_h5($suboutput);
}


sub write_ccs {
  my $bi = shift @_;
  my $ccsoutput = shift @_;
  $bi->get_blasr_ccs_to_genome_alignments_from_h5($ccsoutput);
}
