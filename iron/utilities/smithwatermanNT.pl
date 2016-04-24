#!/usr/bin/perl -w
use strict;
use List::Util qw(max);

### smithwaterman.pl ###
# A very basic Smith-Waterman alignment
# Input: two sequences
# Output: Print a local alignment with the starting nucleotide of each sequence before each sequence
# Modifies: STDOUT

our $maxgap = 10;
our $gapopen = -5;
our $gapextend = -5; 
our $match = 10;
our $mismatch = -5;

#### print_matrix ###
# Input: A 2D array
# Output:  Print matrix to STDOUT
# Modifies:  STDOUT
sub print_matrix {
  my $M = shift @_;
  if(scalar(@{$M})==0) {
    return;
  }
  my $mlen = scalar(@{$M});
  my $nlen = scalar(@{$M->[0]});
  for(my $m = 0; $m < $mlen; $m++) {
    my $oline = '';
    for(my $n = 0; $n < $nlen; $n++) {
      $oline = $oline . ' ' . $M->[$m]->[$n];
    }
    print $oline . "\n";
  }
}

#### diag_score ####
# Fetch the diagnal H value from the current coordinate
# Input:  Matrix H and current coordiante row i and column j
# Output: The Hi-1,j-1 value
# Modifies: None
sub diag_score {
  my $H = shift @_;
  my $i = shift @_;
  my $j = shift @_;
  if($i-1 < 0 or $j-1 < 0) {  return 0; }
  return $H->[$i-1]->[$j-1];
}

#### match_score ####
# Return the score given the current characters
# Input:  Two characters c1 c2
# Output: The score
# Modifies: None
sub match_score {
  my $c1 = shift @_;
  my $c2 = shift @_;
  if($c1 eq $c2) {
    return $match;
  }
  return $mismatch;
}

#### row_scores ####
# Return the scores for the gap going up the row
# Input:  Score matrix H and current position i j
# Output: an array of scores
# Modifies: None
sub row_scores {
  my $H = shift @_;
  my $i = shift @_;
  my $j = shift @_;
  my @oscores;
  if($i==0) {
    push @oscores, 0;
    return @oscores;
  }
  my $bottom = 0;
  if($i-$maxgap > 0) { $bottom = $i-$maxgap; }
  for(my $m = $bottom; $m < $i; $m++) {
    my $k= $i-$m; #distance
    push @oscores, $H->[$i-$k][$j]+$gapopen+($k-1)*$gapextend;
  }
  return @oscores;
}

#### col_scores ####
# Return the scores for the gap going across the columnb
# Input:  Score matrix H and current position i j
# Output: an array of scores
# Modifies: None
sub col_scores {
  my $H = shift @_;
  my $i = shift @_;
  my $j = shift @_;
  my @oscores;
  if($j==0) {
    push @oscores, 0;
    return @oscores;
  }
  my $bottom = 0;
  if($j-$maxgap > 0) { $bottom = $j-$maxgap; }
  for(my $m = $bottom; $m< $j; $m++) {
    my $l=$j-$m; #distance
    push @oscores, ($H->[$i]->[$j-$l]+$gapopen+($l-1)*$gapextend);
  }
  return @oscores;
}

#### score_matrix ###
# Make the H scoring matrix for the alginment
# Input: Two sequences
# Output:  H a matrix with with the scores computed
# Modifies:  STDOUT
sub score_matrix {
  my $s1 = shift @_;
  my $s2 = shift @_;
  #my $H = [[0 for x in range(0,len(s1))] for x in range(0,len(s2))] #initialize alignment matrix
  my @H;
  for(my $i = 0; $i < scalar(@{$s1}); $i++) {
    my @temp;
    $H[$i] = \@temp;
    for(my $j = 0; $j < scalar(@{$s2}); $j++) {
      $H[$i]->[$j] = 0;
    }
  }
  my $mlen = scalar(@H);
  my $nlen = scalar(@{$H[0]});
  for(my $m = 0; $m < $mlen; $m++) {
    for(my $n = 0; $n < $nlen; $n++) { 
      print "$s1->[$m]\t$s2->[$n]\n";
      $H[$m]->[$n] = max(diag_score(\@H,$m,$n) + match_score($s1->[$m],$s2->[$n]),max(row_scores(\@H,$m,$n)),max(col_scores(\@H,$m,$n)),0);
      #print_alignment_matrix(H,s1,s2)
      #print ''
    }
  }
  return \@H;
}

#### matrix_max #####
# return the coordinate of the max value
# Input: takes a matrix H
# Output: list [i,j] with the best coordiante
sub matrix_max {
  my $H = shift @_;
  my $mlen = scalar(@{$H});
  my $nlen = scalar(@{$H->[0]});
  my @best = (0,0);
  my $bestval = 0;
  for(my $m=0; $m < $mlen; $m++) {
    for(my $n=0; $n < $nlen; $n++) {
      if($H->[$m]->[$n] > $bestval) { 
        @best = ($m , $n);
        $bestval = $H->[$m]->[$n]
      }
    }
  }
  return @best;
}

#### next_coord ###
# Print the next coordinate to go to and its score
# Input: H scoring matrix, and current coordinate i j
# Output:  the next score and coordiante inext jnext
# Modifies:  None
sub next_coord {
  my $H = shift @_;
  my $j = shift @_;
  my $i = shift @_;
  my $rowval = 0;
  if($i-1 >= 0) {
    $rowval = $H->[$i-1]->[$j];
  }
  my $colval = 0;
  if($j-1 >= 0) {
    $colval = $H->[$i]->[$j-1];
  }
  my $diagval = 0;
  if($i-1 >=0 and $j-1 >= 0) {
    $diagval = $H->[$i-1]->[$j-1]
  }
  my $iret = $i-1;
  my $jret = $j-1;
  if($diagval >= $rowval && $diagval >= $colval) {
    return ($diagval,$iret,$jret);
  }
  if($rowval >= $colval) {
    return ($rowval,$iret,$j);
  }
  return ($colval,$i,$jret);
}


#### get_local_alignment ###
# Print the local alignment given the scoring matrix
# Input: H scoring matrix, sequences s1 and s2
# Output:  A best local alignment between the two sequences, returns the max alignment score, and two strings that are the alignment lines, and the start indecies for the two sequences, and a descriptor of how s2 differs from s1
# Modifies:  none
sub get_local_alignment {
  my $H = shift @_;
  my $s1 = shift @_;
  my $s2 = shift @_;
  my $mlen = scalar(@{$H});
  my $nlen = scalar(@{$H->[0]});
  my ($i,$j) = matrix_max($H);
  my $currentscore = $H->[$i]->[$j];
  my $maxscore = $currentscore;
  my @s1o;
  my @s2o;
  my $a1 = '';
  my $a2 = '';
  my ($isave,$jsave) = (0,0);
  while($currentscore > 0 && $i >= 0 && $j >= 0) {
    my ($isave, $jsave) = ($i,$j);
    my ($inext,$jnext);
    ($currentscore,$inext,$jnext) = next_coord($H,$i,$j);
    if($inext==$i) { # skip one on s2
      unshift @s1o,$s1->[$i];
      unshift @s2o,'-';
    }
    elsif($jnext==$j) { #skip one on s1
      unshift @s1o, '-';
      unshift @s2o,$s2->[$j];
    } else {
      unshift @s1o,$s1->[$i];
      unshift @s2o,$s2->[$j];
    }
    ($i,$j) = ($inext,$jnext);
  }
  my $s1start = $jsave+1;
  my $s2start = $isave+1;
  $a1 = join('',@s1o);
  $a2 = join('',@s2o);
  return ($maxscore, $a1, $a2,$s1start,$s2start);
}

my $s1 = shift @ARGV;
my $s2 = shift @ARGV;

#M = [[0 for x in range(0,len(s1))] for x in range(0,len(s2))] #initialize alignment matrix
#my @M;
#for(my $i = 0; $i < length($s1); $i++) {
#  my @temp;
#  $M[$i] = \@temp;
#  for(my $j = 0; $j < length($s2); $j++) {
#    $M[$i]->[$j] = 0;
#  }
#}

my @s1 = split(//,$s1);
my @s2 = split(//,$s2);
my $H = score_matrix(\@s1,\@s2);

#print_alignment_matrix(H,s1,s2)
my ($maxscore,$s1align,$s2align,$s1coord,$s2coord) = get_local_alignment($H,\@s1,\@s2);
print $maxscore . "\t" . $s1coord . "\t" . $s1align . "\t" . $s2coord . "\t" . $s2align . "\n";
