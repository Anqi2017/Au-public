#!/usr/bin/perl -w
use strict;
##############
# make a fake fastq file from a fasta file
# use STDIN and STDOUT
# Modifies: File IO

  my $head = "";
  my $seq = "";
  while(my $line = <STDIN>) {
    chomp($line);
    if($line=~/^>(.*)/) { 
      #output the previous
      if($head ne "") {
        print "@"."$head\n$seq\n";
        print "+"."$head\n";
        my $qual = '';
        for(my $i = 0; $i < length($seq); $i++) { $qual .= 'I'; }
        print "$qual\n";
      }
      $head = $1;
      $seq = "";
    } else {
      $seq .= $line
    }
  }
  if($head ne "" and $seq ne "") {
    print "@"."$head\n$seq\n";
    print "+"."$head\n";
    my $qual = '';
    for(my $i = 0; $i < length($seq); $i++) { $qual .= 'I'; }
    print "$qual\n";
  }
