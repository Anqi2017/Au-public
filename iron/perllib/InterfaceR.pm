package InterfaceR;

use strict;

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
  if(-d '/tmp/weirathe') { }
  else { mkdir '/tmp/weirathe'; }
  my $tempdir = '/tmp/weirathe/InterfaceR.'.int(rand()*1000000000);
  mkdir $tempdir;
  if(-d $tempdir) {} else { die "couldn't make temp dir $tempdir\n"; }
  $self->{'tempdir'} = $tempdir;
  return $self;
}

#Split a continuous exposure variable to maximize KM chi
#Otherwise known as the Eharmony Method
#Pre: hash of subjects, each with a surivial time, an event binary (censor), and an exposure (continuous)
#     'time', 'event', and 'exposure' are the keys
#     (optional) a mingroupsize
#Post: An array of low subjects, an array of high subjects, 
#      the range of the cutoff with values for each one
#      a plot showing KM results at various cutoffs
sub MaximumKaplanMeierTimeEventExposure {
  my $self = shift @_;
  my $data = shift @_;
  my $mingroupsize = 1;
  if(scalar(@_) > 0) { $mingroupsize = shift @_; }
  my @pats = sort {$data->{$a}->{'exposure'}<=>$data->{$b}->{'exposure'}} keys %{$data};
  if(scalar(@pats) < 2) { die "too few patients to do maximum KM\n"; }
  my $besti = 0;
  my $bestlow;
  my $besthigh;
  my $bestchi = 0;
  my $bestpval = 0;
  my %all;
  for(my $i = $mingroupsize; $i <= scalar(@pats)-$mingroupsize; $i++) {
    my @km;
    my @low;
    my $loweventcnt = 0;
    for(my $j = 0; $j < $i; $j++) {
      push @low, $pats[$j];
      my %temp;
      $temp{'time'} = $data->{$pats[$j]}->{'time'};
      $temp{'event'} = $data->{$pats[$j]}->{'event'};
      if($temp{'event'}==1) { $loweventcnt++; }
      $temp{'exposure'} = 0;
      push @km, \%temp;
    }
    my @high;
    my $higheventcnt = 0;
    for(my $k = $i; $k < scalar(@pats); $k++) {
      push @high, $pats[$k];
      my %temp;
      $temp{'time'} = $data->{$pats[$k]}->{'time'};
      $temp{'event'} = $data->{$pats[$k]}->{'event'};
      if($temp{'event'}==1) { $higheventcnt++; }
      $temp{'exposure'} = 1;
      push @km, \%temp;
    }
    my $res = $self->KaplanMeierTimeEventExposure(\@km,'pval');
    my %t;
    $all{$i} = \%t;
    $all{$i}->{'chi'} = $res->{'chi'};
    $all{$i}->{'pval'} = $res->{'pval'};
    $all{$i}->{'low'} = scalar(@low);
    $all{$i}->{'high'} = scalar(@high);
    $all{$i}->{'lowevent'} = $loweventcnt;
    $all{$i}->{'highevent'} = $higheventcnt;
    $all{$i}->{'avg_exposure'} = ($data->{$pats[$i-1]}->{'exposure'}+$data->{$pats[$i]}->{'exposure'})/2;
    if($res->{'chi'} > $bestchi) {
      $bestchi = $res->{'chi'};
      $bestpval = $res->{'pval'};
      $besti = $i;
      $bestlow = \@low;
      $besthigh = \@high;
    }
  }
  open(OF,">".$self->{'tempdir'}."/outdata.txt") or die;
  print OF "cutoff\tchi\tavg_exposure\tlowevent\thighevent\n";
  foreach my $cut (sort {$a<=>$b} keys %all) {
    print OF $cut . "\t" . $all{$cut}->{'chi'} . "\t" . $all{$cut}->{'avg_exposure'}."\t".$all{$cut}->{'lowevent'}."\t".$all{$cut}->{'highevent'}."\n";
  }
  close OF;

  open(OF,">".$self->{'tempdir'}."/outplot.R") or die;
  print OF 'd<-read.table("'.$self->{'tempdir'}.'/outdata.txt",header=TRUE)
maxchi = max(d$chi)
minchi = min(d$chi)
maxevent = max(c(max(d$lowevent),max(d$highevent)))
minevent = min(c(min(d$lowevent),min(d$highevent)))
v = (d$avg_exposure-min(d$avg_exposure))*(1/(max(d$avg_exposure)-min(d$avg_exposure)))
v = v*(maxchi-minchi)+minchi
le = (d$lowevent-minevent)/(maxevent-minevent)
le = le*(maxchi-minchi)+minchi
he = (d$highevent-minevent)/(maxevent-minevent)
he = he*(maxchi-minchi)+minchi
png("'.$self->{'tempdir'}."/maxKM.png".'")
plot(cbind(d$cutoff,d$chi),type="n",xlab="low group size",ylab="chi test statistic")
lines(cbind(d$cutoff,d$chi),lwd=3)
lines(cbind(d$cutoff,v),col="blue",lwd=3)
lines(cbind(d$cutoff,le),col="red",lwd=3,lty=2)
lines(cbind(d$cutoff,he),col="yellow",lwd=3,lty=2)
abline(v='.$besti.',col="green",lwd=3)
dev.off()';
  close OF;
  my $cmd = 'Rscript '.$self->{'tempdir'}.'/outplot.R';
  chomp(my @val3s = `$cmd`);

  my %out;
  $out{'chi'} = $bestchi;
  $out{'pval'} = $bestpval;
  $out{'index'} = $besti;
  $out{'low'} = $bestlow;
  $out{'high'} = $besthigh;
  $out{'avg_exposure'} = ($data->{$pats[$besti-1]}->{'exposure'}+$data->{$pats[$besti]}->{'exposure'})/2;
  $out{'all_data'} = \%all;
  $out{'plot_location'} = $self->{'tempdir'}."/maxKM.png";
  return \%out;
}

#Pre: Data for each patient that is time, 
#     whether or not an event has occured, and
#     whether or not the patient has been exposed
#     This can for example by surivial time, whether or not a death or relapse has occured, and whether or not a patient has high expression of a certain gene
#     This data is fed to this function by giving a 
#     per patient or observation array
#     that has 'time', 'event', and 'exposure' fields
#     time is a float, event is binary (1 or 0), and exposure is binary (1 or 0)
#     .. I think event is also called a censor
#     Optional input:  'pval'
#     makes it only output pvalue and chi into the output
sub KaplanMeierTimeEventExposure {
  my $self = shift;
  my $td = $self->{'tempdir'};
  my $pats = shift;
  my $extra = '';
  if(scalar(@_) > 0) { $extra = shift @_; }
  open(OF,">$td/data.txt") or die "could not write $td/data.txt\n";
  print OF "time\tevent\texposure\n";
  foreach my $pat (@{$pats}) {
    print OF $pat->{'time'} . "\t" . $pat->{'event'} . "\t" . $pat->{'exposure'} . "\n";
  }
  close OF;

  open(OF,">$td/KM.R") or die "could not write $td/KM.R\n";

  print OF '#load the KM functions
library("splines")
library("survival")
d<-read.table("'."$td/data.txt".'",header=TRUE)
mfit<-survfit(formula = with(d, Surv(time,event==1)) ~ exposure, data = d)
mdif<-survdiff(Surv(time,event==1) ~ exposure, data = d)';

if($extra ne 'pval') { print OF '
png("'."$td/KM.png".'")
plot(mfit)
dev.off()'; }

print OF  '
print(mdif)';  

  close OF;
  # now we have the data loaded up and read to go lets execute it
  my @vals;
  #chomp(my @vals = `Rscript $td/KM.R`);
  open(OS,"Rscript $td/KM.R |") or die "couldn't execute $td/KM.R.\n1. is R module loaded? 2. can you write there?\n";
  while(my $line = <OS>) {
    chomp($line);
    push @vals, $line;
  }
  close OS;
  my $chi;
  my $pval;
  foreach my $line (@vals) {
    if($line=~/Chisq=\s+(\S+)\s+.*\s+p=\s+(\S+)/) {
      $chi = $1;
      $pval = $2;
    }
  }
  if($extra ne 'pval') {
    open(OF,">$td/report.txt") or die;
    foreach my $line (@vals) {
      print OF "$line\n";  
    }
    close OF;
  }
  my %res;
  $res{'chi'} = $chi;
  $res{'pval'} = $pval;
  $res{'script'} = "$td/KM.R";
  $res{'data_location'} = "$td/data.txt";
  $res{'plot_location'} = "$td/KM.png";
  $res{'report'} = "$td/report.txt";
  return \%res;
}


sub clean {
  my $self = shift @_;
  my $td= $self->{'tempdir'};
  `rm -r $td`;
}

1;
