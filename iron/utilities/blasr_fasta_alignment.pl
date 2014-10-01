#!/usr/bin/perl -w
use strict;
use SequenceBasics qw(rc);
use ErrorCodec qw(rc_code);

# takes to fastas, a subsequence fasta (reads of inserts), and a ccs read in a fasta
if(scalar @ARGV != 2) { die; }
my $subname = shift @ARGV;
my $ccsname = shift @ARGV;
#query is subread
#target is ccs
open(INF,"blasr $subname $ccsname -bestn 100000 -m 0 2>/dev/null |") or die;
my $score = '';
while(my $line = <INF>) {
  chomp($line);
  if($line=~/^\s*Score:\s*(\S+)/) { 
    $score = $1; 
  }
  if($line=~/^\s*Query:\s*(\S+)\/\d+_\d+$/) {
    my $query = $1;
    chomp($line = <INF>);
    my $target = '';
    if($line=~/^\s*Target:\s*(\S+)$/) {
      $target = $1;
    }
    <INF>; #model
    <INF>; #raw score
    <INF>; #map qv
    chomp($line = <INF>); #query strand
    if($line=~/1/) { die; } # unexpected strand direction
    my $strand = '+';
    chomp($line = <INF>); #target strand
    if($line=~/1/) { $strand = '-'; }
    my ($qstart,$qstop,$qlen);
    chomp($line = <INF>); #query range
    if($line=~/^\s*QueryRange:\s*(\d+)\s*->\s*(\d+)\s*of\s*(\d+)$/) {
      ($qstart,$qstop,$qlen) = ($1, $2, $3);
    } else { die; }
    my ($tstart,$tstop,$tlen);
    chomp($line = <INF>); #target range
    if($line=~/\s*TargetRange:\s*(\d+)\s*->\s*(\d+)\s*of\s*(\d+)$/) {
      ($tstart,$tstop,$tlen) = ($1, $2, $3);
    } else { die; }
    my @seqs;
    my $stayin = 1;
    my $tseq = '';
    my $qseq = '';
    while($stayin ==  1) {
      chomp($line=<INF>);
      if(eof(INF) || $line=~/nMatch/) { $stayin = 0; }
      else {
        if($line=~/^\s*$/ && scalar(@seqs) == 2) { #dump lines; 
          $qseq .= $seqs[0];
          $tseq .= $seqs[1];
          #print "$seqs[0]\t$seqs[1]\n";
          @seqs = ();
        } elsif($line=~/^\s*\d+\s*(\S+)/) {
          push @seqs, $1;
        }
      }
    }
    if(scalar(@seqs) == 2) {
      $qseq.=$seqs[0];
      $tseq.=$seqs[1];
    }
    if($strand eq '-') {
      $qseq = rc($qseq);
      $tseq = rc($tseq);
      #$qstart = $qlen - $qstop+1;
      #$tstart = $tlen - $tstop+1;
    } #make the ccs the one oriented correctly
    $qstart++;
    $tstart++; #make index 1 not index 0
    my $ec = new ErrorCodec($tseq,$qseq);
    my $code = $ec->get_error_code();
    if($strand eq '-') {
     $code = rc_code($code);
    }
    print "$target\t$query\t$score\t$strand\t$code\t$tstart\t$qstart\t$tlen\t$qlen\n";
  }
}
