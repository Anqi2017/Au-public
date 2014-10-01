#!/usr/bin/perl -w
use strict;
##### pacbio_raw_to_ccs.pl #####
# Take a .bax.h5 or .bas.h5 and output the <basename>.fasta and <basename>.fastq
# Right now is set with a minimum of 2 passes.
# Input: .bax.h5 or .bas.h5 file name, output base name, accuracy (integer between zero and 100), chemistry (i.e. P4-C2 or unknown)
# Output: <basename>.fasta and <basename>.fastq
# Modifies: STDIO, FileIO, /tmp/weirathe, and maybe smrtanalysis virtual machine may be invoked, but I don't think theres collisions in smrtanalysis's temporary folders
#           Is coded right now for 16 threads, and almost certainly needs to be executed on the cluster to get enough virtual memory

if(scalar(@ARGV) != 4) { die "<infile> <out base name> <accuracy i.e. 90> <chem i.e. P4-C2 or unknown>\n"; }
my $infile = shift @ARGV;
my $outbase = shift @ARGV;
my $accuracy = shift @ARGV;
my $chem = shift @ARGV;
my $rand = int(rand()*10000000);
my $tfolder = "/tmp/weirathe/t$rand";
unless(-d "/tmp/weirathe") {
  `mkdir /tmp/weirathe`;
}
unless(-d "$tfolder") {
  `mkdir $tfolder`;
}
print "$rand\n";
my $moviename;
if($infile=~/([^\/\.]+)\.*\d*\.ba.\.h5$/) {
  $moviename = $1;
} else { die "unrecognized filetype: $infile\n"; }
my $chemname = "$tfolder/chemistry_mapping.xml"; 
my $fofnname = "$tfolder/input.fofn"; 
open(OF,">$chemname") or die; 
print OF "<Map>\n"; 
print OF "  <Mapping><Movie>$moviename</Movie>\n"; 
print OF "  <SequencingChemistry>$chem</SequencingChemistry></Mapping>\n"; 
print OF "</Map>\n"; 
close OF; 
open(OF,">$fofnname") or die; 
print OF "$infile\n"; 
close OF; 
my $cmd1 = '. /Shared/Au/jason/Source/smrtanalysis/current/etc/setup.sh && '; 
$cmd1 .= '/Shared/Au/jason/Source/smrtanalysis/current/analysis/bin/ConsensusTools.sh CircularConsensus '; 
$cmd1 .= '--minFullPasses 2 --minPredictedAccuracy '.$accuracy.' '; 
$cmd1 .= '--chemistry '.$chemname.' '; 
$cmd1 .= '--parameters /Users/weirathe/jason/Source/smrtanalysis/current/analysis/etc/algorithm_parameters/2014-03 '; 
$cmd1 .= '--numThreads 16 --fofn '.$fofnname.' '; 
$cmd1 .= '-o '.$tfolder; 
print "$cmd1\n"; 
open(STR,"$cmd1|") or die; 
while(my $ln = <STR>) {
  chomp($ln);
  print "$ln\n";
}
close STR;
`cp $tfolder/*.fasta $outbase.fasta`;
`cp $tfolder/*.fastq $outbase.fastq`;
`rm -r $tfolder`; 
