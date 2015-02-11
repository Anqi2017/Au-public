#!/usr/bin/perl -w 
use strict;

# take a fasta piped to standard input
# take the number of reads to put in each file
# take an output directory
# will create that directory then put the fasta file broken into peices
# named a by integer for each file.

if(scalar(@ARGV) != 2) { die " cat my.fasta | ./explode_fasta.pl 500 myoutdir/\n"; }
my $num = shift @ARGV;
my $outdir = shift @ARGV;
if(-d $outdir) { die "error $outdir already exists.. won't overwrite.\n"; }
`mkdir $outdir`;
my $seq = '';
my $name = '';
my $printed = $num+1; # force starting a new file
my $filenum = 0;
while(my $line = <STDIN>) {
  chomp($line);
  if($line=~/^(>.+)/) {
    my $oldname = $name;
    $name = $1;
    if($seq ne '') {
      if($printed >= $num) {
        $filenum++;
        close OF;
        open(OF,">$outdir/$filenum.fa") or die;
        $printed = 0;
      }
      $printed++;
      print OF "$oldname\n".$seq."\n";
    }
    $seq = '';
  } else {
    $line=~s/\s//g;
    $seq .= $line;
  }
}
if($seq ne '') {
  if($printed > $num) {
    $filenum++;
    close OF;
    open(OF,"$outdir/$filenum.fa") or die;
    $printed = 0;
  }
  $printed++;
  print OF "$name\n".$seq."\n";
}
close OF;
