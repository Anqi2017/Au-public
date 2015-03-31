#!/usr/bin/perl -w
use strict;
##### pacbio_raw_to_ccs.pl #####
# Input: <ccs fasta> <output base name> <min length> <cpu count> <(optional) custom primers fasta>
#
# Example of custom primer fasta:
#>F0
# 5' sequence here
#>R0
# 3' sequence here (but in reverse complement)
# their explanation of this convention "if you have a full-length read and it has the 5' primer in the first N bases, then one should be able to find an exact match of the 3' primer in the last M bases"
#
# Output: <basename>.fasta and <basename>.fastq
# Modifies: STDIO, FileIO, /tmp/weirathe, and maybe smrtanalysis virtual machine may be invoked, but I don't think theres collisions in smrtanalysis's temporary folders
#           Is coded right now for 16 threads, and almost certainly needs to be executed on the cluster to get enough virtual memory

if(scalar(@ARGV) < 4) { die "<ccs fasta> <output base name> <min length> <cpu count> <(optional) custom primer fasta>\n"; }
my $infile = shift @ARGV;
my $outbase = shift @ARGV;
my $minlength = shift @ARGV;
my $cpus = shift @ARGV;
my $customprimers = "";
if(scalar(@ARGV) == 1) {
  $customprimers = shift @ARGV;
  $customprimers = ' -p '.$customprimers.' ';
}
my $rand = int(rand()*10000000);
my $tfolder = "/localscratch/Users/weirathe/t$rand";
#unless(-d "/tmp/weirathe") {
#  `mkdir /tmp/weirathe`;
#}
unless(-d "$tfolder") {
  `mkdir $tfolder`;
}
print "$rand\n";
my $cmd1 = '. /Shared/Au/jason/Source/smrtanalysis2.3.0/current/etc/setup.sh && '; 
$cmd1 .= '/Shared/Au/jason/Source/smrtanalysis2.3.0/current/analysis/bin/pbtranscript.py classify '; 
$cmd1 .= "$infile $tfolder/temp.fasta ";
$cmd1 .= $customprimers; 
$cmd1 .= ' --cpus '.$cpus.' '; 
$cmd1 .= ' --min_seq_len '.$minlength.' '; 
$cmd1 .= ' -d '."$tfolder/classify".' '; 
$cmd1 .= ' --flnc '."$tfolder/IsoSeqFLNC.fasta "; 
$cmd1 .= ' --nfl '."$tfolder/IsoSeqNFL.fasta"; 
print "$cmd1\n"; 
open(STR,"$cmd1|") or die; 
while(my $ln = <STR>) {
  chomp($ln);
  print "$ln\n";
}
close STR;
`cp $tfolder/IsoSeqFLNC.fasta $outbase.IsoSeqFLNC.fasta`;
`cp $tfolder/IsoSeqNFL.fasta $outbase.IsoSeqNFL.fasta`;
`cp $tfolder/temp.classify_summary.txt $outbase.IsoSeqReport.txt`;
`rm -r $tfolder`; 
