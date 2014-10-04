#!/usr/bin/perl -w
use strict;

# Pre:  1.  A fastq file either gzipped or not.
#           A file with a .gz file extension will be treated as gzipped.
#       2 (optional).  The number of reads to use in the analysis (integer).
#                      Or the random fraction of reads to use (looks for a ".").
# Post: Prints to STDOUT two tables.  The first table, for each 
#       base position, it orders all the quality scores for that position
#       and outputs the quality score at the percentile of that read.
#       so 0% will have the minimum quality 50% will have the median quality and
#       100% will have the maximum quality.
#       The second table is listed per quality score.
#       For each quality score, it a list of cutoffs (number of bases with that quality score or lower).
#       and it lists the percentage of the reads that would be thrown 
#       away if you throw away anything with that quality score or lower for
#       that number of bases.
# Modifies: STDOUT

my $z = 0;
my @positions;
my @qualities;
if(scalar(@ARGV) != 1 && scalar(@ARGV) != 2) { die "./fastq_quality_inspector.pl <fastq file (optionally .gz gzipped)> <optional number of reads (default 10000) as an integer, alternatively a random fraction of reads (i.e. 0.001) when using sorted reads>\n"; }
my $fname = shift @ARGV;
my $opt = '';
if(scalar(@ARGV) > 0) {
  $opt = shift @ARGV;
}
if($fname=~/\.gz$/) { 
  open(INF,"zcat $fname |") or die "failed 'zcat $fname |'\n";
} else {
  open(INF,"$fname") or die;
}
my $stopnum = 10000;
if($opt ne '' && !($opt=~/\./)) { $stopnum = $opt; }
while(my $l1 = <INF>) {
  $z++;
  chomp($l1);
  chomp(my $l2=<INF>);
  chomp(my $l3=<INF>);
  chomp(my $l4=<INF>);
  if($opt=~/\./) { # skip if the optional value says to randomly skip o 
  if($opt < rand()) {
      next; 
    }
  }
  my @chars = split(//,$l4);
  my %qs;
  for(my $i = 0; $i < scalar(@chars); $i++) {
    if(!defined($positions[$i])) {
      my @temp;
      $positions[$i] = \@temp;
    }
    my $onum = ord($chars[$i]);
    if(!exists($qs{$onum})) {
      $qs{$onum} = 0;
    }
    $qs{$onum}++;
    push @{$positions[$i]}, $onum;
  }
  my $totalqs = number_of_bases_equal_or_lower_quality(\%qs);
  push @qualities, $totalqs;
  if($z >= $stopnum  && !($opt=~/\./)) { last; }  
}
close INF;
print_positions(\@positions);
print_qualities(\@qualities);

sub print_qualities {
  my $qualities = shift @_;
  my %qs;
  my $maxcount=0;
  my $maxq = 0;
  my $minq = 9**9**9;
  foreach my $read (@{$qualities}) {
    foreach my $q (keys (%{$read})) {
      if($q < $minq) { $minq = $q; } 
      if($q > $maxq) { $maxq = $q; } 
      if($read->{$q} > $maxcount) { $maxcount = $read->{$q}; }
    }
  }
  my %s;
  for(my $q = $minq; $q <= $maxq; $q++) {
    my %temp;
    $s{$q} = \%temp;
    for(my $i = 1; $i <= $maxcount; $i++) {
      $s{$q}->{$i} = 0;
    }
  }
  # fill in gaps
  foreach my $read (@{$qualities}) {
    for(my $q = $minq; $q <= $maxq; $q++) {
      if(!exists($read->{$q})) {
        for(my $j = $minq; $j < $q; $j++) {
          if(exists($read->{$j})) {
            $read->{$q} = $read->{$j};
          }
        }
      }
      if(!exists($read->{$q})) {
       $read->{$q} = 0;
      }
    }
  }
  foreach my $read (@{$qualities}) {
    foreach my $q (keys %{$read}) {
      my $cnt = $read->{$q};
      for(my $basecount = 1; $basecount <= $cnt; $basecount++) {
        $s{$q}->{$basecount}++;
      }
    }
  }
  my $readcount = scalar(@{$qualities});
  print "Quality\tPercent of reads with X number of bases of Quality or lower.\n";
  print "-------\t-----------------------------------------------------------\n";
  foreach my $q (sort{$a<=>$b} keys %s) {
    print "$q\t";
    my $ostring = '';
    foreach my $c (sort {$a<=>$b} keys %{$s{$q}}) {
      my $val = sprintf("%.3g",$s{$q}->{$c}*100/$readcount);
      $ostring .= "$c:$val\%,";
    }
    chop($ostring);
    print "$ostring\n";
  }
  print "-------\t-----------------------------------------------------------\n";
  print "Quality\tPercent of reads with X number of bases of Quality or lower.\n";
  print "-------\t-----------------------------------------------------------\n";
  print "used $readcount reads\n";
}

sub deep_copy {
  my $ref = shift @_;
  my %out;
  foreach my $k (keys %{$ref}) {
    $out{$k} = $ref->{$k};
  }
  return \%out;
}

sub add_to_qualities {
  my $qualities = shift @_;
  my $vals = shift @_;
  my %qnums;
  foreach my $q (keys %{$qualities}) {
    $qnums{$q} = 0;
  }
  foreach my $q (keys %{$vals}) {
    $qnums{$q} = 0;
  }
  my @qs = sort {$a<=>$b} keys %qnums;
  my $min = $qs[0];
  my $max = $qs[scalar(@qs)-1];
  @qs = ();
  for(my $i = $min; $i <= $max; $i++) {
    push @qs, $i;
  }
  for(my $i = 0; $i < scalar(@qs); $i++) {
    # fill in missing main
    if(!exists($qualities->{$qs[$i]})) { 
      foreach my $q (@qs) {
        if(exists($qualities->{$q}) && $q < $qs[$i]) {
          $qualities->{$qs[$i]} = deep_copy($qualities->{$q});
        }
      }
    }
    if(!exists($qualities->{$qs[$i]})) { 
      my %temp;
      $qualities->{$qs[$i]} = \%temp;
    }
    # fill in missing alts
    if(!exists($vals->{$qs[$i]})) {
      foreach my $q (@qs) {
        if(exists($vals->{$q}) && $q < $qs[$i]) {
          $vals->{$qs[$i]} = $vals->{$q};
        }
      }
    }
    if(!exists($vals->{$qs[$i]})) {
      $vals->{$qs[$i]} = 0;
    }
    # now this quality exists in both sets for sure
    # for this quality want to say for each base count cutoff how many reads are included
    my $cnt = $vals->{$qs[$i]};
    #print "q: $qs[$i] $cnt\n";
    if(!exists($qualities->{$qs[$i]}->{$cnt})) {
      $qualities->{$qs[$i]}->{$cnt} = 0;
    }
    $qualities->{$qs[$i]}->{$cnt}+=1;
  }
  return;
}

sub number_of_bases_equal_or_lower_quality {
  my $ref = shift @_;
  my @qs = sort {$a<=>$b} keys %{$ref};
  my %out;
  for(my $i = 0; $i < scalar(@qs); $i++) {
    my $sum = 0;
    for(my $j = 0; $j <= $i; $j++) {
      $sum += $ref->{$qs[$j]};
    }
    $out{$qs[$i]} = $sum;
  }
  return \%out;
}

sub print_positions {
  my $positions = shift @_;
  my @total;
  my @perc = (1,5,10,25,50,75,90,95,99);
  print "Pos\t0,";
  foreach my $j (@perc) {
    print "$j,";
  }
  print "100%\n";
  print "-----\t---------------------------------\n";
  for(my $i = 0; $i < scalar(@{$positions}); $i++) {
    my $pos = $i+1;
    print "$pos\t";
    my @vals = sort {$a<=>$b} @{$positions[$i]};
    push @total, @vals;
    print "$vals[0],";
    foreach my $j (@perc) {
      my $v = percentile(\@vals,$j);
      print "$v,";
    }
    print "$vals[scalar(@{$positions[$i]})-1]\n";
  }
  print "-----\t---------------------------------\n";
  print "Pos\t0,";
  foreach my $j (@perc) {
    print "$j,";
  }
  print "100%\n";
  print "-----\t---------------------------------\n";
  @total = sort {$a<=>$b} @total;
  print "Total\t";
  print "$total[0],";
  foreach my $j (@perc) {
    my $v = percentile(\@total,$j);
    print "$v,";
  }
  print "$total[scalar(@total)-1]\n";
  print "\n";  
}

sub percentile {
  my $ref = shift @_;
  my $num = shift @_;
  return $ref->[scalar(@{$ref})*$num/100];
}
