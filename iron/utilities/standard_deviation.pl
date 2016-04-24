#!/usr/bin/perl -w
use strict;

#Pre: User pipes a stream of numbers into stdin
#Post: Prints average and standard deviation for entry
#      so last printed number is average and standard deviation for everything

my $sum = 0;
my $cnt = 0;
my $sumofsquares = 0;
while(my $line = <STDIN>) {
  chomp($line);
  unless($line=~/^[\d+\-\.]+$/) { die "non numeric: $line\n"; }
  $sum += $line;
  $cnt++;
  my $mean = $sum/$cnt;
  $sumofsquares += ($line-$mean)*($line-$mean);
  my $s = 'NA';
  if($cnt > 1) {
    $s = sqrt($sumofsquares/($cnt-1));
  }
  print "$mean\t$s\n";
}
