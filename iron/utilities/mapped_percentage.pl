#!/usr/bin/perl -w
use strict;

# What percentage of reads mapped to a file?
# Pre: Stream a samfile in from STDIN, and take a list of fastq files
# Post: Output to standard output read counts in, read counts out, and percent mapped.
# Modifies: Standard output

if(scalar(@ARGV) < 1) { die; }
my $incnt = 0;
foreach my $file (@ARGV) {
  open(INF,$file) or die;
  while(my $line = <INF>) {
    if($line=~/^@/) { $incnt++; }
  }
  close INF;
}
my $outcnt = 0;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^@/) { }
  else { $outcnt++; }
}

#print "Input: $incnt\n";
#print "Output: $outcnt\n";
my $p = $outcnt/$incnt*100;
my $ln = sprintf("%.2f",$p);
print "$ln\n";
