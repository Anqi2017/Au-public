#!/usr/bin/perl -w
use strict;
##############
# make a fasta file from a fastq file
# if no input file, use STDIN and STDOUT
# Input:  fastq file name (optional), fasta file name (optional)
# Output:  fasta file
# Modifies: File IO

#if(scalar(@ARGV) == 0) { die "<input filename (optional)> <output filename (optional)>\n"; }
# usage: <input filename (optional)> <output filename (optional)>
while(my $line1 = <STDIN>) {
    chomp($line1);
    my $name;
    if($line1=~/^@(.+)/) { $name = $1; } else { die "bad header1 format\n"; }
    if($name=~/[\t]/) { die "unsupported tab character in sequence name\n"; }
    chomp(my $line2 = <STDIN>);
    chomp(my $line3 = <STDIN>);
    if($line3=~/^\+/) { } else { die "bad header2 format\n"; }
    chomp(my $line4 = <STDIN>);
    print "$name\t$line2\t$line3\t$line4\n";
}
