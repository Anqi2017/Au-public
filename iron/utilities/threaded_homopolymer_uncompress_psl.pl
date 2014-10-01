#!/usr/bin/perl -w
use strict;
use HomopolymerCompression;
use SequenceBasics qw(read_fasta_into_hash rc);
use EnhancedAligner;
require threads;

require threads::shared;

if(scalar(@ARGV) != 5) { die "<input .psl (hp compressed)> <output .psl> <uncompressed genome (fasta)> <uncompressed long reads (fasta)> <thread count>\n"; }

my $inpsl = shift @ARGV;
my $outpsl = shift @ARGV;
my $original_genome = shift @ARGV;
my $original_query = shift @ARGV;
my $thread_count = shift @ARGV;

my $qseqs_1 = read_fasta_into_hash($original_query);
## deep copy
my %intemp : shared;
foreach my $t (keys %{$qseqs_1}) {
  $intemp{$t} = $qseqs_1->{$t};
}
my $qseqs : shared = \%intemp;

## find the chromosomes you need to do first
open(INF,$inpsl) or die;
my %chroms;
while(my $line = <INF>) {
  chomp($line);
  my @fields = split(/\t/,$line);
  $chroms{$fields[13]}=  1;
}
close INF;

open(OF,">$outpsl") or die;
foreach my $chrom (sort {$b cmp $a} keys %chroms) {
  my $psl_entries : shared = get_chromosome_from_psl($chrom,$inpsl);
  print "read psl entires\n";
  my $thc = get_chromosome_from_genome($chrom,$original_genome);
  print "got chromosome $chrom compressed and uncompressed.\n";
  my $uncompressed_chromosome : shared = $thc->get_uncompressed_nucleotides();
  my $tSize = length($thc->get_uncompressed_nucleotides());
  my $numthreads = $thread_count;
  my $conversion_map : shared = $thc->get_compressed_to_uncompressed_conversion_map();
  print "got conversion map for $chrom\n";
  my @threads;
  print $chrom . "\n";
  for(my $i = 0; $i < $numthreads; $i++) { 
    my $thr = threads->create( \&do_psls, $psl_entries,$qseqs,$chrom,$uncompressed_chromosome,$conversion_map,$tSize,$i,$numthreads);
    push @threads, $thr;
    #my $line = do_psls($psl_entries,$qseqs,$chrom,$thc,$tSize,$i,$numthreads);
    #print $line;
  }
  my @results;
  foreach my $thr (@threads) {
    push @results, $thr->join();
  }
  foreach my $result (@results) { print OF $result; }
}
close OF;

sub local_convert {
  my $map = shift @_;
  my $coord1 = shift @_;
  my $coord2 = shift @_;
  $coord1--;
  $coord2--;
  my ($i1,$temp1) = @{$map->[$coord1]};
  my ($temp2,$i2) = @{$map->[$coord2]};
  return ($i1+1,$i2+1);
}

sub do_psls {
  print "tdoing thread\n";
  #my $thc;
  my $psl_entries = shift @_;
  my $qseqs=shift @_;
  my $chrom = shift @_;
  my $uncompressed_chromosome = shift @_; # bring in the chromsome
  my $conversion_map = shift @_;
  my $tSize = shift @_; # get the size too
  my $thread_index = shift @_;
  my $total_threads = shift @_;
  my $oline = '';
  my $z = 0;
  foreach my $psl (@{$psl_entries}) {
   $z++;
   if($z%$total_threads==$thread_index) { 
    my @fields = split(/\t/,$psl);
    my $qname = $fields[9];
    my $tname = $fields[10];
    my $strand = $fields[8];
    my $qhc = new HomopolymerCompression();
    if(!exists($qseqs->{$fields[9]})) { die "missing $fields[9]\n"; }
    if($strand eq '-') {
      $qhc->load_uncompressed_nucleotides(rc($qseqs->{$fields[9]}));
    } else {
      $qhc->load_uncompressed_nucleotides($qseqs->{$fields[9]});
    }
    my $qSize = length($qhc->get_uncompressed_nucleotides());
    $qhc->create_map_compressed_to_uncompressed();
    my ($sizes,$qstarts,$tstarts) = ($fields[18],$fields[19],$fields[20]);
    my @sizes = split(/\s*,\s*/,$sizes);
    my @qstarts = split(/\s*,\s*/,$qstarts);
    my @tstarts = split(/\s*,\s*/,$tstarts);
    my $name = $fields[13];
    my $windowsize = 10;
    ## adding one to each of these because of that weird dumb half open index zero format of psl 
    #my ($iprevstart,$iprevend) = $thc->convert_compressed_range_to_uncompressed_range_coordinate($tstarts[0]+1,$tstarts[0]+$sizes[0]-1);
    my ($iprevstart,$iprevend) = local_convert($conversion_map,$tstarts[0]+1,$tstarts[0]+$sizes[0]-1);
    my ($jprevstart,$jprevend) = $qhc->convert_compressed_range_to_uncompressed_range_coordinate($qstarts[0]+1,$qstarts[0]+$sizes[0]-1);
    my @query_coords;
    my @target_coords;
    for(my $i = 0; $i < scalar(@sizes); $i++) {
      my $i1 = $tstarts[$i]+1; #again this is for that weird half index zero stuff
      my $i2 = $tstarts[$i]+$sizes[$i]-1;
      my $j1 = $qstarts[$i]+1; #half index zero again thanks
      my $j2 = $qstarts[$i]+$sizes[$i]-1;
      #my ($ui1,$ui2) = $thc->convert_compressed_range_to_uncompressed_range_coordinate($i1,$i2);
      my ($ui1, $ui2) = local_convert($conversion_map,$i1,$i2);
      if($ui1 > $iprevend + $windowsize) { #output old set
        my @temp = ($iprevstart, $iprevend);
        push @target_coords, \@temp;
        $iprevstart = $ui1;
      } 
      $iprevend = $ui2;
      my ($uj1,$uj2) = $qhc->convert_compressed_range_to_uncompressed_range_coordinate($j1,$j2);
      if($uj1 > $jprevend + $windowsize) { #output old set
        my @temp = ($jprevstart, $jprevend);
        push @query_coords, \@temp;
        $jprevstart = $uj1;
      } 
      $jprevend = $uj2;
    }
    my @temp1 = ($iprevstart, $iprevend);
    push @target_coords, \@temp1;
    my @temp2 = ($jprevstart, $jprevend);
    push @query_coords, \@temp2;
    #### now we have target and query coordinates;
    my @target_map;
    foreach my $coord (@target_coords) {
      #my $seq = $thc->get_uncompressed_range($coord->[0],$coord->[1]);
      my $seq = substr($uncompressed_chromosome,$coord->[0]-1,$coord->[1]-$coord->[0]+1);
      my @chars = split(//,$seq);
      my $h = 0;
      for(my $i = $coord->[0]; $i <= $coord->[1]; $i++) {
        my @temp = ($chars[$h],$i);
        push @target_map, \@temp;
        $h++;
      }
    }
    my @query_map;
    foreach my $coord (@query_coords) {
      my $seq = $qhc->get_uncompressed_range($coord->[0],$coord->[1]);
      my @chars = split(//,$seq);
      my $h = 0;
      for(my $i = $coord->[0]; $i <= $coord->[1]; $i++) {
        my @temp = ($chars[$h],$i);
        push @query_map, \@temp;
        $h++;
      }
    }
    #target map and query map are ready for alignment
    my $tar = '';
    foreach my $c (@target_map) { $tar .= $c->[0]; }
    my $que = '';
    foreach my $c (@query_map) { $que .= $c->[0]; }
    ######## This is where to realign.
    #if($strand eq '+') {
      #print "$tar\n$que\n\n\n";
      my $ea = new EnhancedAligner($tar,$que);
      my ($ta,$qa)=$ea->homopolymer_needleman_wunsch();
      #print "$ta\n$qa\n\n";
      #print seqlen($ta) . "\t" . length($tar) . "\n";
      my @tacoords;
      my @ta = split(//,$ta);
      my $k = 0;
      my $mismatches = 0;
      foreach(my $i = 0; $i < scalar(@ta); $i++) {
        if($ta[$i] ne '-') { 
          push @tacoords,$target_map[$k]->[1];
          #print $target_map[$k]->[0] . "\t" . $target_map[$k]->[1] . "\n";
          $k++; 
        } else {
          push @tacoords, '-';
        }
      }
      my @qacoords;
      my @qa = split(//,$qa);
      my $m = 0;
      foreach(my $i = 0; $i < scalar(@ta); $i++) {
        if($qa[$i] ne '-') { 
          push @qacoords,$query_map[$m]->[1];
          #print $query_map[$m]->[1] . "\n";
          $m++; 
        } else {
          push @qacoords, '-';
        }
      }
      my $matches = 0;
      for(my $i = 0; $i < scalar(@qa); $i++) {
        if($qa[$i] ne '-' && $ta[$i] ne '-') { 
          if($ta[$i] ne $qa[$i]) { $mismatches++; } 
          else { $matches++; }
        }
      }
      my @both;
      for(my $i = 0; $i < scalar(@qacoords); $i++) {
        if($tacoords[$i] ne '-' && $qacoords[$i] ne '-') {
          my @temp = ($tacoords[$i],$qacoords[$i]);
          push @both,\@temp;
        }
      }
      my $pt = $both[0]->[0];
      my $pq = $both[0]->[1];
      @qstarts = (); #saving over that original i guess
      my @qends;
      @tstarts = ();
      my @tends;
      #print "start $pt\t$pq\n";
      push @tstarts,$pt;
      push @qstarts,$pq;
      for(my $i = 1; $i < scalar(@both); $i++) {
        my @coord = @{$both[$i]};
        if($coord[0] > $pt+1 || $coord[1] > $pq + 1 && $pt > 0) {
          #we've jumped more than one and its greater than zero so flush the buffer
          #print "end $pt\t$pq\n";
          push @tends,$pt;
          push @qends,$pq;
          #print "start $coord[0]\t$coord[1]\n";
          push @tstarts,$coord[0];
          push @qstarts,$coord[1];
        }
        $pt = $coord[0];
        $pq = $coord[1];
      }
      push @tends, $both[scalar(@both)-1]->[0];
      push @qends, $both[scalar(@both)-1]->[1];
      my $blocksizes='';
      $qstarts=''; #saving over the original i guess
      $tstarts='';
      my $qStart = $qstarts[0]-1;
      my $qEnd = $qends[scalar(@tends)-1];
      my $tStart = $tstarts[0]-1;
      my $tEnd = $tends[scalar(@tends)-1];
      for(my $i = 0; $i < scalar(@qstarts); $i++) {
        $blocksizes.=($tends[$i]-$tstarts[$i]+1) . ',';
        $qstarts.=($qstarts[$i]-1).','; #that dumb half index zero thing
        $tstarts.=($tstarts[$i]-1).',';
      }
      my $qNumInsert = 0;
      my $qBaseInsert = 0;
      my $tNumInsert = 0;
      my $tBaseInsert = 0;
      for(my $i = 1; $i < scalar(@qstarts); $i++) {
        if($qstarts[$i]-$qends[$i-1] > 1) { 
          $qNumInsert++;
          $qBaseInsert+=$qstarts[$i]-$qends[$i-1];
        }
        if($tstarts[$i]-$tends[$i-1] > 1) { 
          $tNumInsert++;
          $tBaseInsert+=$tstarts[$i]-$tends[$i-1];
        }
      }
      my $toline = ''; 
      $toline.= "$matches\t$mismatches\t0\t0\t$qNumInsert\t$qBaseInsert\t";
      $toline.= "$tNumInsert\t$tBaseInsert\t$strand\t$qname\t$qSize\t$qStart\t$qEnd\t";
      $toline.= "$chrom\t$tSize\t$tStart\t$tEnd\t$blocksizes\t$qstarts\t$tstarts\n";
      #print $toline;
      $oline.=$toline;
    #} doing for either strand now        
    ######## Now we know the alignments, we map coordinates to each base
   }
  }
  return $oline;
}

sub seqlen {
  my $seq = shift @_;
  $seq=~s/-//g;
  return length($seq);
}

sub get_chromosome_from_psl {
  my $chrom = shift @_;
  my $inpsl = shift @_;
  my @lines : shared;
  open(INF,$inpsl) or die;
  while(my $line = <INF>) {
    chomp($line);
    my @fields = split(/\t/,$line);
    if($fields[13] eq $chrom) {
      push @lines,$line;
    }
  }
  return \@lines;
}

sub get_chromosome_from_genome {
  my $chrom = shift @_;
  my $original_genome = shift @_;
  open(INF,$original_genome) or die;
  my $head='';
  my $seq='';
  my %seqs;
  while(my $line = <INF>) {
    chomp($line);
    if($line=~/^>\s*(\S+)/) { #try the first nonwhitespace element as the title
      if($seq ne '' && $head eq $chrom) {
        #print "compressing $head\n";
        my $hc = new HomopolymerCompression();
        $hc->load_uncompressed_nucleotides($seq);
        #print "creating maps for $head\n";
        $hc->create_map_compressed_to_uncompressed();
        close INF;
        return $hc;
      }
      $head = $1;
      $seq = '';
    } else {
      $seq .= $line;
    }
  }
  if($seq ne '' && $head eq $chrom) {
    #print "compressing $head\n";
    my $hc = new HomopolymerCompression();
    $hc->load_uncompressed_nucleotides($seq);
    #print "creating maps for $head\n";
    $hc->create_map_compressed_to_uncompressed();
    close INF;
    return $hc;
  }
  die "$chrom was not found in the reference genome\n";
  close INF;
}
