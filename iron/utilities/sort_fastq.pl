#!/usr/bin/perl -w
use strict;
#Pre: A fastq file streamed in STDIN
#Post: A fastq file ordered by the sequence name streamed to STDOUT

open(STREAM,"| fastq_to_tsv.pl | sort -k 1,1 | tsv_to_fastq.pl > /dev/stdout") or die;
while(my $line = <STDIN>) {
  print STREAM $line;
}
