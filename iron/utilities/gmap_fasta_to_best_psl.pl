#!/usr/bin/perl -w
use strict;
use SequenceBasics qw (best_psl);
##### gmap_fasta_to_best_psl.pl #############
# Launch a very fast (optionally multiple threaded) long read 
#        aligner that returns a sorted psl file with only the 
#        best alignment per query.
# Input: A fasta file to map to the genome
#        A filename for where to write psl file
#        A path to GMAP index dataset (a directory)
#        (optional) A number of threads (like -t)
# Output: A psl file.
# Requirements: samtools binaries AND MISC scripts 
#               must be in the path
#               gmap binaries must be in the path
# Modifies: The output file for the sorted bam,
#           not sure if sorting makes a temp anywhere


my $t = 1;
if(scalar(@ARGV) < 3) { die "<input fasta> <output (*.psl)> <gmap index path (directory)> <number of threads (optional)>\n"; }
my $fname = shift @ARGV;
my $fout = shift @ARGV;
my $gmapindexpath = shift @ARGV;
if($gmapindexpath=~/^(.*)\/$/) { $gmapindexpath=$1; } 
my $gmapindexname;
if($gmapindexpath=~/([^\/]+)$/) {
  $gmapindexname=$1;
} else { die "trouble geting index name\n"; }
if(scalar(@ARGV) == 1) {
  $t = shift @ARGV;
}
#my $gmap_cmd = 'gmap -D '.$gmapindexpath.' -f 1 -d '.$gmapindexname.' -t '.$t.' '.$fname.' 2>/dev/null';
my $gmap_cmd = 'gmap -D '.$gmapindexpath.' -f 1 -d '.$gmapindexname.' -t '.$t.' '.$fname;
my $rnum = int(10000000*rand());
my $tdir = "/tmp/weirathe.$rnum";
`mkdir $tdir`;
#my $cmd = $gmap_cmd.' | psl2sam.pl | samtools view -Sb -t '.$genomeindexpath.' - | samtools sort - '.$fout;
my $cmd = $gmap_cmd." > $tdir/all.psl";
print "$cmd\n";
`$cmd`;
best_psl("$tdir/all.psl",$fout);
`rm -r $tdir/`;
