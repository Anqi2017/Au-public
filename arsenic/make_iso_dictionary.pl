#!/usr/bin/perl -w
# -w provides warnings

#This script makes a reference table to associate isoforms with gene
#names and subscrips.  It takes a refFlat file as an argument.

use strict; 


my $refflat_file = shift @ARGV;


open (INF, $refflat_file) or die;

my %genes;
my %isos;

while (my $refflatline = <INF>)
{
	chomp ($refflatline); 
	if ($refflatline =~ /^#/) {next;} #skips header

	my @fields = split(/\t/,$refflatline);
	my $gene = $fields[0];
	my $iso = $fields[1];

	if (!exists($genes{$gene})) #if new gene, make a reference 
	{
		my %temp;
		$genes{$gene} = \%temp;
	}
	
	$genes{$gene} -> {$iso} = 1; #gene references iso
}


foreach my $gene(keys %genes)
{
	my @isos = keys(%{$genes{$gene}});
	my $i = 0;
	my $size = scalar(@isos);

	for ($size)
	{
		if ($size == 1)
		{
			print($gene, "\t", $isos[0], "\t", $gene, "\n");
		}
		
		elsif ($size > 1)
		{
			for ($i = 0; $i < $size; $i++)
			{
				print ($gene, "\t", $isos[$i], "\t", $gene, ".", $i, "\n");
			}		
		}
	}
}
