#!/usr/bin/perl -w 
use strict;

# take a fasta piped to standard input
# write a fasta to stdout if length of sequence is longer than input length

my $seq = '';
my $name = '';
my $inlen;
if(scalar (@ARGV) != 1) { die "cat seq.fasta | ./filter_fasta_by_length.pl > seq2.fasta\n"; }
$inlen = shift @ARGV;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^>([^\t]+)/) {
    my $oldname = $name;
    $name = $1;
    if($seq ne '' && length($seq) >= $inlen) {
      print ">$oldname\n$seq\n";
    }
    $seq = '';
  } else {
    $line=~s/\s//g;
    $seq .= $line;
  }
}
if($seq ne '' && length($seq) >= $inlen) {
  print ">$name\n$seq\n";
}
