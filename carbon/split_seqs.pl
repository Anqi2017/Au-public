#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $fastq_flag = '';
my $max_length;
my $tsv_flag = '';
GetOptions("length=i" => \$max_length,
           "fastq"    => \$fastq_flag, 
           "tsv"      => \$tsv_flag
) or die("Error in command line arguments\n");


my $line;
my @a; 
my $length;
my @b; 
my $quot;
my $rem;

while($line = <>){
  chomp $line;
  @a = split(/\t/,$line);
  $length = length($a[1]);
  @b = ();
  if($length > $max_length){
    $quot = int( $length / $max_length);
    $rem  = $length % $max_length; 
    for(my $i=0; $i < $quot; $i++){
      if( $fastq_flag ){
        if( $tsv_flag){
          push(@b, "$a[0]_$i\t".substr($a[1], $i*$max_length, $max_length)."\t$a[2]\t".substr($a[3],$i*$max_length, $max_length)); 
        }else{
          push(@b, "$a[0]_$i\n".substr($a[1], $i*$max_length, $max_length)."\n$a[2]\n".substr($a[3],$i*$max_length, $max_length)); 
        }
      }else{
        if( $tsv_flag){
          push(@b, "$a[0]_$i\t".substr($a[1], $i*$max_length, $max_length)); 
        }else{
          push(@b, ">$a[0]_$i\n".substr($a[1], $i*$max_length, $max_length)); 
        }
      }
    }
    if( $rem > 0 ){
      if( $fastq_flag){
        if( $tsv_flag){
         push(@b, "$a[0]_".$quot."\t".substr($a[1], $max_length*$quot, $rem)."\t$a[2]\t".substr($a[3],$quot*$max_length, $rem));
        }else{
         push(@b, "$a[0]_".$quot."\n".substr($a[1], $max_length*$quot, $rem)."\n$a[2]\n".substr($a[3],$quot*$max_length, $rem));
        }
      }else{
        if( $tsv_flag){
         push(@b, "$a[0]_".$quot."\t".substr($a[1], $max_length*$quot, $rem));
        }else{
         push(@b, ">$a[0]_".$quot."\n".substr($a[1], $max_length*$quot, $rem));
        }
      }
    }
    print join("\n",@b), "\n";
  }else{
    if( $fastq_flag ){
      if( $tsv_flag){
        print "$a[0]\t$a[1]\t$a[2]\t$a[3]\n";
      }else{
        print "$a[0]\n$a[1]\n$a[2]\n$a[3]\n";
      }
    }else{
      if( $tsv_flag){
        print "$a[0]\t$a[1]\n";
      }else{
        print ">$a[0]\n$a[1]\n";
      }
    }
  }
}
