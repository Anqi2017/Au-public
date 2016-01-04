#!/usr/bin/perl -w
# -w provides warnings

#This script replaces the isoforms in a table with names defined in an exernal
#table of format "gene"     \t      "isoform"      \t      "new name".
#It takes the table to be changed as the first argument and the names table
#as the second argument.

use strict; 

my $iso_file = shift @ARGV;

my $ref_file = shift @ARGV;

my %nameof;


open (INF, $ref_file) or die;
while (my $refline = <INF>)
{
	chomp ($refline);
	
	my @fields = split(/\t/,$refline);
	#my $gene = $fields[0];
	my $iso = $fields[1];
	my $name = $fields[2];

	$nameof{$iso} = $name;
}

close ($ref_file);

open (INF2, $iso_file) or die;
my $header = <INF2>;
print ($header);
while (my $line = <INF2>)
{
	chomp ($line);
	my @fields = split(/\t/,$line);
	my $iso = $fields[0];
	my $size = scalar(@fields);

	if (exists $nameof{$iso})
	{
		print ($nameof{$iso});
	}

	else 
	{
		print ($iso);
	}

	for (my $i = 1; $i < $size; $i++)
	{
		print("\t", $fields[$i]);
	}
	print ("\n");
}

close ($iso_file);






