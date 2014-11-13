#!/usr/bin/perl -w
# -w provides warnings

use strict; 

my $cmd = "ls /Shared/Au/jason/Code/HoustonJack_20140830/HUMANCUFFREFSEQ/*/genes*";

chomp(my @src_files = `$cmd`);
#@ is an array

my %genes;
# % designates a hash

my %samples;

foreach my $gene_file(@src_files)
{
	#print "$gene_file\n";
	my $sample;
	if ($gene_file =~/\/BL(\d+)\./) {
	#forward slashes surround regular expressions
	#backslash is necessary to have a forward slash inside it
	#\d is 1 character 0-9, \d+ catches all such characters
	#. would mean 1 of anything without \
	#() designate what will be $1

		$sample = $1;
	}
	else { die; }

	$samples{$sample} = 1;	

	open (INF,$gene_file) or die;
	my $header = <INF>;
	while (my $line = <INF>) {
		chomp ($line);
		my @fields = split(/\t/,$line);
		if ($fields[6] =~/_/) {
			next;
		}
		my $gene_name = $fields[0];
		my $fpkm = $fields[9];
		if (!exists($genes{$gene_name})) {
			my %temp;
			$genes{$gene_name}=\%temp;
			# \ to reference
		}
               
		if(!exists($genes{$gene_name}->{$sample})) {
			$genes{$gene_name}->{$sample} = 0;
		}
		$genes{$gene_name}->{$sample}+=$fpkm;
	}

}

my @ordered_samples = sort {$a<=>$b} keys %samples;

print "#Sample";

foreach my $sample (@ordered_samples)
{
	print "\t$sample";
}

print "\n";

foreach my $gene_name(keys %genes)
{
	print "$gene_name\t";
	my $buffer = '';

	foreach my $sample (@ordered_samples) {
		$buffer .= $genes{$gene_name}->{$sample}."\t";
		#.= means string concatenation
	}
	chop ($buffer);
	#removes last character

	print "$buffer\n";
}
