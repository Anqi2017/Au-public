#!/usr/bin/perl -w
use strict;
#### launcher for picard to make it more like sam tools
my $picdir = '/Shared/Au/jason/Source/picard-tools-1.112';
chomp(my @progs = `ls $picdir/*.jar`);
my %progs;
foreach my $prog (@progs) {
  my $name;
  if($prog =~/([^\/]+).jar/) { $name = $1; } else { die "funkyname\n"; }
  $progs{$name} = $prog;
}
if(scalar(@ARGV) == 0) {
  foreach my $name (sort {$a cmp $b} keys %progs) {
    print $name ."\n";
  }
  print "try one of the above programs with its required arguments.";
} else {
  my $inname = shift @ARGV;
  if(scalar(@ARGV) == 0) {
    if(exists($progs{$inname})) {
      chomp(my @lines = `java -Xmx2g -jar $progs{$inname} --help`);
      foreach my $line (@lines) { print $line . "\n"; }
    } else {
      print "unrecognized picard program\n";
    }
  } else {
    my $cmd = "java -Xmx2g -jar $progs{$inname}";
    foreach my $arg (@ARGV) {
      $cmd .= " $arg";
    }
    open(INF,"$cmd"."|") or die "error launching $cmd\n";
    while(my $line = <INF>) {
      print $line;
    }
    close INF;    
  }
}
