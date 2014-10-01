#!/usr/bin/perl -w
use strict;
use threads;
use SequenceBasics qw(break_fasta best_psl);

##### blat_fasta_to_best_psl.pl #############
# Launch a fast (optionally multiple threaded) long read 
#        aligner that returns a sorted bam file. Yes.
# Input: A fasta file to map to the genome
#        A filename for where to write the psl file
#        A path to BLAT database
#        (optional) A number of threads (like -t)
# Output: A psl file.
# Requirements: BLAT binary (blat) must be in path 
# Modifies: The output file for the psl


my $t = 1;
if(scalar(@ARGV) < 3) { die "<input (query) fasta> <output (*.psl)> <blat database (*.fa or *.nib or *.2bid)> <number of threads (optional)>\n"; }
my $rnum = int(rand()*10000000);
my $tdir = "/tmp/weirathe.$rnum";
`mkdir $tdir`;
print "made\n$tdir\ntemporary directory\n\n";

my $fname = shift @ARGV;
my $fout = shift @ARGV;
if($fout=~/^(.*)\.bam$/) { $fout = $1; }
my $blatindexpath = shift @ARGV;
if(scalar(@ARGV) == 1) {
  $t = shift @ARGV;
}
break_fasta($fname,"$tdir/part.fa",$t);
my @bthreads;
for(my $i = 1; $i <= $t; $i++) {
  my $th = threads->create(\&execute_blat,$blatindexpath,"$tdir/part.fa.$i","$tdir/part.psl.$i");
  push @bthreads, $th;
}
foreach my $th (@bthreads) { $th->join(); }
for(my $i = 1; $i <= $t; $i++) {
  `cat $tdir/part.psl.$i >> $tdir/all.psl`;  
}
best_psl("$tdir/all.psl","$fout");
`rm -r $tdir`;

sub execute_blat {
  my $blatindexpath = shift @_;
  my $fasta = shift @_;
  my $output = shift @_;
  my $blat_cmd = 'blat '.$blatindexpath.' '."$fasta".' '."$output";
  print "$blat_cmd\n";
  `$blat_cmd`;
}
