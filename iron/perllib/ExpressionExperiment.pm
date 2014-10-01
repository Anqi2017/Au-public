package ExpressionExperiment;

use strict;

use Scalar::Util qw(looks_like_number);

use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw();
%EXPORT_TAGS = ();

sub new {
  my $self = {};
  my $class = shift @_;
  bless $self;
  my @temp;
  $self->{'data'} = \@temp;
  my %temp;
  $self->{'probes'} = \%temp;
  return $self;
}

sub get_name {
  my $self = shift @_;
  return $self->{'name'};
}

sub is_set {
  my $self = shift @_;
  if(scalar(@_) != 1) { die "is_set(probe)\n"; }
  my $probe = shift @_;
  my $ind = $self->{'probes'}->{$probe};
  my $val = $self->{'data'}->[$ind];
  if($val eq '') { return 0; }
  return 1;  
}

sub get_expression {
  my $self = shift @_;
  if(scalar(@_) != 1 ) {  die "get_expression(probe)\n"; }
  my $probe = shift @_;
  my $ind = $self->{'probes'}->{$probe};
  my $val = $self->{'data'}->[$ind];
  if($val eq '') { return; }
  return $val;
}

#Post: output an array of probe names in the order they were given
sub get_probe_list {
  my $self = shift @_;
  return sort {$self->{'probes'}->{$a}<=>$self->{'probes'}->{$b}} keys %{$self->{'probes'}};
}

sub set_name {
  my $self = shift @_;
  if(scalar(@_) != 1) { die "set_name(inputname)\n"; }
  $self->{'name'} = shift @_;
}

sub add_measurement {
  my $self = shift @_;
  if(scalar(@_) != 2) { die "add_measurement(probe,value)\n"; }
  my $probe = shift @_;
  my $value = shift @_;
  if(exists($self->{'probes'}->{$probe})) {
    die "already have a measurement for $probe, it takes unqiuely named probes for each experiment\n";
  }
  if(!looks_like_number($value)) {
    print STDERR "$probe has a value $value which should be a number but doesn't look like one.\n";
    $value = '';
  }
  my $i = scalar(@{$self->{'data'}});
  $self->{'probes'}->{$probe} = $i;
  push @{$self->{'data'}}, $value;
}

1;
