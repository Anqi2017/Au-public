#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

my $fastq_flag;
my $max_length;
GetOptions("fastq=i"  => \$fastq_flag, 
           "length=i" => \$max_length
) or die("Error in command line arguments\n");


my $line;
my @a; 
my $length;
my @b; 
my $quot;
my $rem;

my @seq;
while($line = <>){
  chomp $line;
  @a = split(/\t/,$line);

  if($fastq_flag){
    die("Not yet implemented\n");
  }else{
    @seq = split("", $a[1]);
    $length = length($a[1]);
    #print STDERR "ON: $a[0], max length: $max_length, length: $length\n";
    if($length > $max_length){
      $quot = int( $length / $max_length);
      $rem  = $length % $max_length; 
      for(my $i=0; $i < $quot; $i++){
        #print STDERR "SPLITTING: $a[0], iter $i\n";
        #print STDERR substr($a[1], $i*$max_length, $max_length), "\n";
        push(@b, "$a[0]_$i\t".substr($a[1], $i*$max_length, $max_length)); 
      }
      if( $rem > 0 ){
        push(@b, "$a[0]_".$quot."\t".substr($a[1], $max_length*$quot, $rem));
      }
      print join("\n",@b), "\n";
    }else{
      print "$a[0]\t$a[1]\n";
    }

  }
}
