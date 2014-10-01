package ExpressionExperimentSet;

use strict;
use Digest::MD5 qw(md5_hex);

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
  $self->{'experiments'} = \%temp;
  my %temp2;
  $self->{'probes'}=\%temp2;
  return $self;
}

sub get_experiment {
  my $self = shift @_;
  if(scalar(@_) != 1) { die "get_experiment(experiment_name)\n"; }
  my $name = shift @_;
  return $self->{'data'}->[$self->{'experiments'}->{$name}];
}
#Pre: a probe name
#Post: a hash of expression values keyed by experiment name
#      only defined values are included
sub get_expression_values_by_probe {
  my $self = shift @_;
  if(scalar(@_) != 1) { die "get_expression_values_by_probe(probe)\n"; }
  my $probe = shift @_;
  my @experiment_list = $self->get_experiment_name_list();
  my %res;
  foreach my $exname (@experiment_list) {
    my $ex = $self->get_experiment($exname);
    if($ex->is_set($probe))  {  $res{$ex->get_name()} = $ex->get_expression($probe); }
  }
  return \%res;
}

#Post: return an array of experiment names in the order they were input
sub get_experiment_name_list {
  my $self = shift @_;
  return sort {$self->{'experiments'}->{$a}<=>$self->{'experiments'}->{$b}} keys %{$self->{'experiments'}};
}

sub get_probe_list {
  my $self = shift @_;
  my @experiments = $self->get_experiment_name_list();
  my %md5vals;
  my @probeset;
  foreach my $experimentname (@experiments) {
    my $ex = $self->get_experiment($experimentname);
    @probeset = $ex->get_probe_list();
    my $md5val = '';
    foreach my $probe (@probeset) {
      $md5val = md5_hex($md5val . $probe);
    }
    $md5vals{$md5val} = 1;
  }
  if(scalar(keys %md5vals) > 1) { die "Different probe list between experiments\nPlease access them within each experiment if they are right.. it could just be they are in a different order\n"; }
  return @probeset;
}

#faster version of get_probe_list when you're sure they are all the same
sub get_probe_list_greedy {
  my $self = shift @_;
  my @experiments = $self->get_experiment_name_list();
  my $ex = $self->get_experiment($experiments[0]);
  my @probeset = $ex->get_probe_list();
  return @probeset;
}

sub add_experiment {
  my $self = shift @_;
  if(scalar(@_) != 1) { die "add_experiment(ExpressionExperiment)\n"; }
  my $experiment = shift @_;
  if(exists($self->{'experiments'}->{$experiment->get_name()})) { die "experiment $experiment has already been added\n"; }
  my $i = scalar(@{$self->{'data'}});
  $self->{'experiments'}->{$experiment->get_name()} = $i;
  push @{$self->{'data'}}, $experiment;
}

# Pre:  An expression table where the columns have experiment name labels at the top
#       and the rows have probe names, and everything else is expression level
#       File is a TSV (tab seperated value)
sub read_expression_table {
  my $self = shift @_;
  if(scalar(@_) != 1) { die "read_expression_table(expressiontable.txt)\n"; }
  my $filename = shift @_;
  open(INF,"$filename") or die "could not open $filename\n";
  chomp(my $heads = <INF>);
  print STDERR "starting reading expression table\n";
  my @heads = split(/\t/,$heads);
  shift @heads;
  foreach my $experiment (@heads) {
    my $ex = ExpressionExperiment->new();
    $ex->set_name($experiment);
    $self->add_experiment($ex);
  }
  my $j = 0;
  while(my $line = <INF>) {
    $j++;
    if($j%1000==0) { print STDERR "$j lines read\n"; }
    chomp($line);
    my @fields = split(/\t/,$line);
    my $probe = shift @fields;
    if(scalar(@fields) != scalar(@heads)) { die "incomplete row error $j\n"; }
    for(my $i = 0; $i < scalar(@fields); $i++) {
      my $experiment = $heads[$i];
      my $expression = $fields[$i];
      $self->get_experiment($experiment)->add_measurement($probe,$expression);
    }
  }
  close INF;
  print STDERR "done reading expression table\n";
}

1;
