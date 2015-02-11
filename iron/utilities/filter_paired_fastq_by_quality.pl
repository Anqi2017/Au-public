#!/usr/bin/perl -w
use strict;
# Pre: 1., 2. Takes two paired end fastqs as inputs can be gzipped if they end
#             in .gz, in this case zcat must be present in the path.
#      3., 4. Takes two output file names.  These will not zip as its written 
#             now.
#      5.  The quality score at which or be low you want to remove sequences
#      6.  The number of bases as or below your quality threshold needing to be
#          present before you remove a sequence from the file.
# Post:       Writes the filtered fastq to the two output files.  
#             Two more files will be created with your outputfilenames+'.bad' 
#             that contain the removed sequences.  It will not overwrite files
#     ----->  If a read failed in one of the of the two mate pairs, both mates
#             are removed.
# Modifies:  Creates files.  Creates those '.bad' files based on your 
#            output file names.  
if(scalar(@ARGV) != 6) { die "filter_fastq_based_on_quality.pl <input1> <input2> <output1> <output2> <low quality value to cut on or below> <low quality cutoff count cut on or above>\n"; }
my $input1 = shift @ARGV;
my $input2 = shift @ARGV;
my $output1 = shift @ARGV;
my $output2 = shift @ARGV;
my $qvalue = shift @ARGV;
my $qcount = shift @ARGV;
my $isgzip = 0;
if($input1=~/\.gz$/) {
  open(INF1,"zcat $input1 |") or die;
  $isgzip = 1;
} else {
  open(INF1,"$input1") or die;
}
if($input2=~/\.gz$/) {
  open(INF2,"zcat $input2 |") or die;  
} else {
  open(INF2,"$input2") or die;
}
if(-e "$output1") { die "$output1 already exists\n"; }
if(-e "$output2") { die "$output1 already exists\n"; }
if(-e "$output1.bad") { die "$output1.bad already exists\n"; }
if(-e "$output2.bad") { die "$output1.bad already exists\n"; }
if($output1=~/\.gz$/ && $output2=~/\.gz$/) {
  open(OF1,"| gzip >$output1") or die;
  open(OF2,"| gzip >$output2") or die;
  open(OF1b,"| gzip >$output1.bad.gz") or die;
  open(OF2b,"| gzip >$output2.bad.gz") or die;
} else {
  open(OF1,">$output1") or die;
  open(OF2,">$output2") or die;
  open(OF1b,">$output1.bad") or die;
  open(OF2b,">$output2.bad") or die;
}
while(my $a1 = <INF1>) {
  chomp($a1);
  chomp(my $a2=<INF1>);
  chomp(my $a3=<INF1>);
  chomp(my $a4=<INF1>);
  chomp(my $b1=<INF2>);
  chomp(my $b2=<INF2>);
  chomp(my $b3=<INF2>);
  chomp(my $b4=<INF2>);
  my @chars1 = split(//,$a4);
  my $lowcount1 = 0;
  foreach my $char (@chars1) {
    if(ord($char) <= $qvalue) {
      $lowcount1++;
    }
  }
  my @chars2 = split(//,$b4);
  my $lowcount2 = 0;
  foreach my $char (@chars2) {
    if(ord($char) <= $qvalue) {
      $lowcount2++;
    }
  }
  if($qcount <= $lowcount1 || $qcount <= $lowcount2) {
    print OF1b "$a1\n$a2\n$a3\n$a4\n";
    print OF2b "$b1\n$b2\n$b3\n$b4\n";
  } else {
    print OF1 "$a1\n$a2\n$a3\n$a4\n";
    print OF2 "$b1\n$b2\n$b3\n$b4\n";
  }
}
close INF1;
close INF2;
