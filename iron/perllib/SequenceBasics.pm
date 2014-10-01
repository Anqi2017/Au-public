package SequenceBasics;

use strict;
require threads;
require threads::shared;

use List::Util qw(shuffle);

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(rc homopolymer_shuffle break_fasta break_fasta_by_max_len best_psl random_dna_sequence read_fasta_into_array read_fasta_into_hash get_sequence_from_coordinates);
%EXPORT_TAGS = ();



##### random_dna_sequence ######
# Make a random sequence of a length
# Input: A length
# Output: A string random ATCG nucletide
# Modifies: None
sub random_dna_sequence {
  my $length = shift @_;
  my $o = '';
  for(my $i = 0; $i < $length; $i++) {
    my $rnum = rand();
    if($rnum < 0.25) {
      $o .= 'A';
    } elsif($rnum < 0.5) {
      $o .= 'T';
    } elsif($rnum < 0.75) {
      $o .= 'C';
    } else {
      $o .= 'G';
    }
  }
  return $o;
}


##### homopolymer_shuffle ######
# Randomize a sequence by shuffeling homopolymers
# Input: Sequence string
# Output: Sequence string
# Modifies: None
sub homopolymer_shuffle {
  my $seq = shift @_;
  my @chars = split(//,$seq);
  my $prev = '';
  my @homos;
  foreach my $char (@chars) {
    if($char ne $prev) {
      push @homos,$char;
    } else {
      $homos[scalar(@homos)-1].=$char;
    }
    $prev = $char;
  }
  my @rand = shuffle(@homos);
  return join('',@rand);  
}


##### rc ######
# Just reverse complement the sequence
# Input: sequence string
# Output: reverse complemented sequence string
# Modifies: none
sub rc { 
  my $seq = shift @_;
  if($seq=~/[uU]/ && ($seq=~/[tT]/)) { die "Mix of Uu Tt in sequence.  I don't know if its RNA or DNA." }
  my $rna = 0;
  if($seq=~/[uU]/) {
    $rna = 1;
  }
  my @chars = split(//,$seq);
  my $o = '';
  foreach my $c (@chars) {
    if($c eq 'A' && $rna == 0) { $o = 'T' . $o; }
    elsif($c eq 'A' && $rna == 1) { $o = 'U' . $o; }
    elsif($c eq 'C') { $o = 'G' . $o; }
    elsif($c eq 'G') { $o = 'C' . $o; }
    elsif($c eq 'T') { $o = 'A' . $o; }
    elsif($c eq 'U') { $o = 'A' . $o; }
    elsif($c eq 'a' && $rna == 0) { $o = 't' . $o; }
    elsif($c eq 'a' && $rna == 1) { $o = 'u' . $o; }
    elsif($c eq 'c') { $o = 'g' . $o; }
    elsif($c eq 'g') { $o = 'c' . $o; }
    elsif($c eq 't') { $o = 'a' . $o; }
    elsif($c eq 'u') { $o = 'a' . $o; }
    else { $o = $c . $o; } #just use this odd character as is.
  }
  return $o;
}

#### best_psl ####
# Take a psl format file
# and make a new psl file with the best hit for each read
sub best_psl {
  my $infile = shift @_;
  my $outfile = shift @_;
  open(INF,$infile) or die "couldn't open psl\n";
  my %best;
  my %bestnum;
  while(my $line = <INF>) {
    chomp($line);
    my @fields = split(/\t/,$line);
    if(scalar(@fields) == 21) { #if its not part of a header
      if($fields[0]=~/^\d+$/) { #if its not part of a header
        if(!exists($bestnum{$fields[9]})) {
          $best{$fields[9]} = $line;
          $bestnum{$fields[9]} = $fields[0];
        } elsif($fields[0]>$bestnum{$fields[9]}) {
          $best{$fields[9]} = $line;
          $bestnum{$fields[9]} = $fields[0];
        }
      }
    } 
  }
  open(OF,">$outfile") or die "couldn't open output psl\n";
  foreach my $m (sort {$a cmp $b} keys %best) {
    print OF "$best{$m}\n";
  }
  close INF;
  close OF;
}

#### break_fasta_by_max_len ###
# input: input fasta, output basename, maxlength
# output: broken files basename.1 basename.2
# modifies: File IO
sub break_fasta_by_max_len {
  my $infasta = shift @_;
  my $outbase = shift @_;
  my $maxlen = shift @_;
  my $entries = read_fasta_into_array($infasta);
  my $i = 0;
  my $j = 0;
  foreach my $entry (@{$entries}) {
    if($i%$maxlen==0) {
      $j++;
      close OF;
      open(OF,">$outbase.$j") or die;
    }
    print OF $entry->{'header'} . "\n" . $entry->{'seq'} . "\n";
    $i++;
  }
}

#### break_fasta ####
# input: input fasta, output basename, number of files to make
# output: broken files filebase.1 filebase.2 filebase.2
# modifies: File IO
sub break_fasta {
  my $infile = shift @_;
  my $outbase = shift @_;
  my $pieces = shift @_;
  my @handles;
  for(my $i = 1; $i <= $pieces; $i++) {
    my $thandle;
    open($thandle,">$outbase.$i") or die;
    push @handles, $thandle;
  } 
  open(INF,$infile) or die;
  my $buffer = '';
  my $i = 0;
  while(my $line = <INF>) {
    chomp($line);
    if($line=~/^>/) {
      if($buffer ne '') { #if something is there
        my $t = $i%$pieces;
        my $f = $handles[$t];
        print $f $buffer;
        $i++;
        $buffer = '';
      }
    }
    $buffer .= "$line\n";
  }
  if($buffer ne '') { 
    my $t = $i%$pieces;
    my $f = $handles[$t];
    print $f $buffer; 
  }
  foreach my $handle (@handles) { close $handle; }
}

#### read_fasta_into_array
#  Input: fasta file
#  Output: Array of hashes with a 'seq' and 'header' key for each entry
sub read_fasta_into_array {
  my $filename = shift @_;
  open(INF,$filename) or die;
  my $head = '';
  my $seq = '';
  my @seqs;
  while(my $line = <INF>) {
    chomp($line);
    if($line=~/^>(.*)/) {
      if($seq ne '') {
        my %temp;
        $temp{'seq'} = $seq;
        $temp{'header'} = $head;
        push @seqs, \%temp;
      }
      $head = $1;
      $seq = '';
    } else { $seq .= $line; }
  }
  if($seq ne '') {
    my %temp;
    $temp{'seq'} = $seq;
    $temp{'header'} = $head;
    push @seqs, \%temp;
  }
  close INF;
  return \@seqs;
}

#### read_fasta_into_hash ###
# Input: fasta file
# Output: Hash with headers as keys and sequence as the value
# Warning: Overwrites identical header entires.  
sub read_fasta_into_hash {
  my $filename = shift @_;
  open(INF,$filename) or die "failed to open $filename for reading\n";
  my $head = '';
  my $seq = '';
  my %seqs;
  while(my $line = <INF>) {
    chomp($line);
    if($line=~/^>(.*)/) {
      if($seq ne '') {
        ### ignore repeats #### if(exists($seqs{$head})) { die "fasta fead file failure:  repeat entry in fasta file\n"; }
        $seqs{$head}=$seq;
      }
      $head = $1;
      $seq = '';
    } else { $seq .= $line; }
  }
  if($seq ne '') {
    $seqs{$head}=$seq;
  }
  close INF;
  return \%seqs;
}

#Pre: genome_hash, $chromosome, 1 indexed start, 1 indexed finish (like ucsc or bed)
#Post: sequence
sub get_sequence_from_coordinates {
  my $fasta = shift @_;
  my $chr = shift @_;
  my $start = shift @_;
  my $stop = shift @_;
  my $seq = '';
  if(!exists($fasta->{$chr})) { 
    die "Trying to access chromosome $chr not in fasta reterning empty string\n";
  }
  if(length($fasta->{$chr}) < $stop) {
    die  "Too short coordinate $chr for $stop stop\n";
  }
  $seq = substr($fasta->{$chr},$start-1,$stop-$start+1);
  return $seq;  
}

1;
