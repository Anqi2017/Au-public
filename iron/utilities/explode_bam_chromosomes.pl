#!/usr/bin/perl -w 
use strict;
################
# Break a bam into chromosomes
# Input: bam filename, outputbasename
# Output: outputbasename is the prefix, then .1.bam incremented
# Modifies: FileIO, stderr i think, forks to # of chr, creates a bam index in the same directory as the input bam if one isn't already there

if(scalar(@ARGV) < 2) { die "<inputbam> <outputbasename>\n"; }
my $inbam = shift @ARGV;
my $outname = shift @ARGV;
open(INF,"samtools view $inbam | cut -f 3 |") or die;
my %chrs;
while(my $chr = <INF>) {
  chomp($chr);
  $chrs{$chr} = 1;
}
close INF;
unless(-e "$inbam.bai") {
  `samtools index $inbam`;
}
my $i = 0;
foreach my $chr (keys %chrs) {
  $i++;
  my $pid = fork();
  if($pid) {}
  elsif($pid==0) {
    sleep(rand()*5);
    my $cmd = "samtools view -bh $inbam $chr | samtools sort -o - - > $outname.$i.bam";
    print "$cmd\n";
    `$cmd`;
    exit 0;
  } else { die "fork fail\n"; }
}
1 while wait() >= 0;
print "finished explode\n";
