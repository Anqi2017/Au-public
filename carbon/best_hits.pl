#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $file_type;
my $header = '';
GetOptions( "type=s"   => \$file_type,
            "header" => \$header
);

my $line;
my %hash;
my @a;
my @b; 
while($line =<>){
  chomp $line;
  if( $header and $. == 1){
    next;
  }
  @a = split(/\t/,$line);

  if( $type eq "blast" ){
    if( exists($hash{$a[0]}) ){
      @b = split(/\t/,$hash{$a[0]});
      if( $a[2] > $b[2]){
        $hash{$a[0]} = $line;
      }
    }else{
      $hash{$a[0]} = $line;
    }
  }elsif($type eq "psl" ){

  }elsif($type eq "lastz-tab" ){
    if( exists($hash{$a[6]}) ){
      @b = split(/\t/,$hash{$a[0]});
      if( $a[0] > $b[0] ){
        $hash{$a[6]} = $line;
      }
    }else{
      $hash{$a[6]} = $line;
    }
  }elsif($type eq "blasr-tab" ){
    if( exists($hash{$a[0]}) ){
      @b = split(/\t/,$hash{$a[0]});
      if( $a[2] > $b[2] ){
        $hash{$a[0]} = $line;
      }
    }else{
      $hash{$a[0]} = $line;
    }
  }else{
    die("File type not allowed. Must be one of [blast,psl,lastz-tab,blasr-tab]\n");
  }
}

foreach my $val (values %hash){
  print "$val\n";
}
