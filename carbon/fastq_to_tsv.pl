#!/usr/bin/perl -w 
use strict;

# take a fastq piped to standard input
# write a tsv to standard output

my $header;
my $seq;
my $header2;
my $quality;
while(my $line = <STDIN>) {
  chomp $line;
  if( $. % 4 == 1){
    $header = $line;
  }elsif( $. % 4 == 2){
    $seq = $line;
  }elsif( $. % 4 == 3){
    $header2 = $line;
  }elsif($. % 4 == 0){
    $quality = $line;
    print "$header\t$seq\t$header2\t$quality\n";
    $header='';
    $seq='';
    $header2='';
    $quality='';
  }
}
