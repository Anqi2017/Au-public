#!/usr/bin/perl -w
use strict;
#Pre: A fasta file streamed in STDIN
#Post: A fasta file ordered by the sequence name streamed to STDOUT

open(STREAM,"| fasta_to_tsv.pl | sort -k 1,1 | tsv_to_fasta.pl > /dev/stdout") or die;
while(my $line = <STDIN>) {
  print STREAM $line;
}
