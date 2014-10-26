#!/usr/bin/perl
use warnings;
use strict;

my $line;
my %hash;
my @a;
my @b; 
while($line =<>){
  chomp $line;
  @a = split(/\t/,$line);
  if( exists($hash{$a[0]}) ){
    @b = split(/\t/,$hash{$a[0]});
    if( $a[2] > $b[2]){
      $hash{$a[0]} = $line;
    }
  }else{
    $hash{$a[0]} = $line;
  }
}

foreach my $val (values %hash){
  print "$val\n";
}
