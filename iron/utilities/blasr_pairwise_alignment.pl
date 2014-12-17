#!/usr/bin/perl -w
use strict;
#Pre: Two sequences 
#    s1 should be expected sequence (target)
#    s2 should be observed sequence (query)
#Post: alignent
#  <orientation> its - if we had to flip the target to match the query
#  <n matches> the number of matches
#  <target alignment>
#  <query alignment>
if(scalar @ARGV != 2) { die; }

my $s1 = shift @ARGV;
my $s2 = shift @ARGV;
my $tdir2 = "/tmp/weirathe";
unless(-d "$tdir2") {
  `mkdir $tdir2`;
}
my $rnum = int(rand()*10000000);
# takes to fastas, a subsequence fasta (reads of inserts), and a ccs read in a fasta
my $s1name = "$tdir2/blasr.1.$rnum.fa";
my $s2name = "$tdir2/blasr.2.$rnum.fa";
open(OF1,">$s1name") or die;
print OF1 ">read1\n$s1\n";
close OF1;
open(OF2,">$s2name") or die;
print OF2 ">read2\n$s2\n";
close OF2;
#query is subread
#target is ccs
open(INF,"blasr $s2name $s1name -bestn 100000 -m 2 2>/dev/null |") or die;
my $res = '';
while(my $line = <INF>) {
  chomp($line);
  $res .= $line;
}
close INF;
my $strand;
if($res=~/targetStrand="([^"]*)"/) {
  $strand = $1;
}
print "$strand\t";
if($res=~/<\s*nCorrect value="([^"]*)"/) {
  print $1."\t";
}
if($res=~/<\s*target\s*>(.*)<\s*target\s*\/>/) {
  print $1."\t";
}
if($res=~/<\s*query\s*>(.*)<\s*query\s*\/>/) {
  print $1."\n";
}
`rm $s1name`;
`rm $s2name`;
