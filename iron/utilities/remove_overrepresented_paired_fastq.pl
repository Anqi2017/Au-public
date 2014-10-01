#!/usr/bin/perl -w
use strict;
#use Digest::SHA1 qw(sha1_hex);

###################
# remove overrepresented sequences from a fastq
# Pre:  a fastq file and choose a threshold for percentage that is acceptable, requires samtools
#       can be a gziped fastq
# Post: a fastq file without that sequence
# Modifies: None

my $cnt = 0;
my %seqs;
my %seqs_b;
my $tot = 0;
my $tot_b = 0;
if(scalar(@ARGV) < 6) { die "<infile_1> <infile_2> <outfile_1> <outfile_2> <threshold percent> <length to look>\n"; }
my $filename1 = shift @ARGV;
my $filename2 = shift @ARGV;
my $outfile1 = shift @ARGV;
my $outfile2 = shift @ARGV;
our $threshold = shift @ARGV;
our $look = shift @ARGV;
if($filename1=~/\.gz$/) {
  open(INF1,"zcat $filename1 |") or die;
} else {
  open(INF1,"$filename1") or die;
}
if($filename2=~/\.gz$/) {
  open(INF2,"zcat $filename2 |") or die;
} else {
  open(INF2,"$filename2") or die;
}
while(my $line1 = <INF1>) {
  $tot++;
  $tot_b++;
  chomp($line1);
  chomp(my $line1_b = <INF2>);
  chomp(my $line2 = <INF1>);
  chomp(my $line2_b = <INF2>);
  #my $val = sha1_hex(substr($line2,0,$look));
  my $val = substr($line2,0,$look);
  if(!exists($seqs{$val})) {
    $seqs{$val} = 0;
  }
  $seqs{$val}++;
  #my $val_b = sha1_hex(substr($line2_b,0,$look));
  my $val_b = substr($line2_b,0,$look);
  if(!exists($seqs_b{$val_b})) {
    $seqs_b{$val_b} = 0;
  }
  $seqs_b{$val_b}++;
  chomp(my $line3 = <INF1>);
  chomp(my $line3_b = <INF2>);
  chomp(my $line4 = <INF1>);
  chomp(my $line4_b = <INF2>);
}
close INF1;
close INF2;

my @over;
my @over_b;

my $ntot = get_over(\%seqs,\@over,$tot);
my $ntot_b = get_over(\%seqs_b,\@over_b,$tot_b);

while($ntot < $tot) {
  $tot = $ntot;
  $ntot = get_over(\%seqs,\@over,$tot);
}

while($ntot_b < $tot_b) {
  $tot_b = $ntot_b;
  $ntot_b = get_over(\%seqs_b,\@over_b,$tot_b);
}
my %bad;
my %bad_b;
foreach my $over (@over) {
  $bad{$over} = 1;
}
foreach my $over_b (@over_b) {
  $bad_b{$over_b} = 1;
}

if($filename1=~/\.gz$/) {
  open(INF1,"zcat $filename1 |") or die;
} else {
  open(INF1,"$filename1") or die;
}
if($outfile1=~/\.gz$/) {
  open(OF1,"| gzip > $outfile1") or die;
} else {
  open(OF1,">$outfile1") or die;
}

if($filename2=~/\.gz$/) {
  open(INF2,"zcat $filename2 |") or die;
} else {
  open(INF2,"$filename2") or die;
}
if($outfile2=~/\.gz$/) {
  open(OF2,"| gzip > $outfile2") or die;
} else {
  open(OF2,">$outfile2") or die;
}
if($outfile1=~/^(.*)\.gz$/) {
  open(OFB1,"| gzip > $1.bad.gz") or die;
} else {
  open(OFB1,">$outfile1.bad") or die;
}
if($outfile2=~/^(.*)\.gz$/) { 
  open(OFB2,"| gzip > $1.bad.gz") or die;
} else {
  open(OFB2,">$outfile2.bad") or die;
}
while(my $line1 = <INF1>) {
  chomp($line1);
  chomp(my $line2 = <INF1>);
  chomp(my $line3 = <INF1>);
  chomp(my $line4 = <INF1>);
  chomp(my $line1_b = <INF2>);
  chomp(my $line2_b = <INF2>);
  chomp(my $line3_b = <INF2>);
  chomp(my $line4_b = <INF2>);
  my $set = "$line1\n$line2\n$line3\n$line4";
  my $set_b = "$line1_b\n$line2_b\n$line3_b\n$line4_b";
  #my $val = sha1_hex(substr($line2,0,$look));
  my $val = substr($line2,0,$look);
  #my $val_b = sha1_hex(substr($line2_b,0,$look));
  my $val_b = substr($line2_b,0,$look);
  if(exists($bad{$val})) {
    print OFB1 "$set\n"; 
  } elsif(exists($bad_b{$val_b})) { 
    print OFB2 "$set_b\n"; }
  else { 
    print OF1 "$set\n"; 
    print OF2 "$set_b\n"; 
  }
}
close OF1;
close OFB1;
close INF1;
close OF2;
close OFB2;
close INF2;


sub get_over {
  my $seqref = shift @_;
  my $overref = shift @_;
  my $tot = shift @_;
  foreach my $seq (keys %{$seqref}) {
    my $num = $seqref->{$seq};
    #print "$seq\n$num\t$tot\n\n";
    if($num/$tot > $threshold/100) {
      push @{$overref}, $seq;
      $tot-=$num;
      delete $seqref->{$seq};
    }
  }
  return $tot;
}
