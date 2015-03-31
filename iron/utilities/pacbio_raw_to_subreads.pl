#!/usr/bin/perl -w
use strict;
##### pacbio_raw_to_subreads.pl ########
# Convert a .bas.h5 or .bax.h5 file into a fasta and fastq
# Input: <input file name>, <output file name base>
# Output: <output file name base>.fasta and <output file name base>.fastq
# Modifies:  Temporary files in /tmp/, STDIO, File IO and anything smrtanalysis is doing
#            Probably can only be called on the cluster since the VMs invoked by smrtanalysis won't be able to launch on the headnodes

if(scalar(@ARGV) != 2) { die; }
my $infile = shift @ARGV;
my $basefile = shift @ARGV;

my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);

#unless(-d "/tmp/$username") {
#  `mkdir /tmp/$username`;
#}
my $rand = int(rand()*100000000);
my $path = "/localscratch/Users/$username/t".$rand;
unless(-d "$path") {
  `mkdir $path`;
}
unless(-d "$path/mytemp") {
  `mkdir $path/mytemp`;
}
my $txml = $path."/input.xml";
print "$txml\n";
open(OF,">$txml") or die;
print OF '<?xml version="1.0"?>
<pacbioAnalysisInputs>
  <dataReferences>
    <url ref="run:0000000-0000"><location>'.$infile.'</location></url>
  </dataReferences>
</pacbioAnalysisInputs>' . "\n";
close OF;

my $tsettings = $path."/settings.xml";
print "$tsettings\n";
open(OF,">$tsettings") or die;
print OF '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<smrtpipeSettings>
    <protocol version="2.3.0" id="RS_Subreads.1" editable="false">
    <application>Data Prep</application>
        <param name="name" label="Protocol Name">
            <value>RS_Subreads</value>
            <input type="text"/>
            <rule required="true" message="Protocol name is required"/>
        </param>
        <param name="description">
            <value>Filter subreads based on read length and quality, optionally splitting by barcode. Output FASTA and bas.h5 file of subreads.</value>
            <textarea></textarea>
        </param>
        <param name="version" hidden="true">
            <value>1</value>
            <input type="text"/>
            <rule type="digits" required="true" min="1.0"/>
        </param>
        <param name="state">
            <value>active</value>
            <input value="active" type="radio"/>
            <input value="inactive" type="radio"/>
        </param>
        <param name="control" hidden="true">
            <value></value>
        </param>
        <param name="fetch" hidden="true">
            <value>common/protocols/preprocessing/Fetch.1.xml</value>
        </param>
        <param name="filtering">
            <value>common/protocols/filtering/SFilter.1.xml</value>
            <select multiple="true">
                <import extension="xml" contentType="text/directory">common/protocols/filtering</import>
            </select>
        </param>
        <param name="barcode" editableInJob="true">
            <value>common/protocols/barcode/NoBarcode.1.xml</value>
            <select multiple="false">
                <import extension="xml" contentType="text/directory">common/protocols/barcode</import>
            </select>
        </param>
    </protocol>
    <moduleStage name="fetch" editable="true">
        <module label="Fetch v1" id="P_Fetch" editableInJob="true">
            <description>Sets up inputs</description>
        </module>
    </moduleStage>
    <moduleStage name="filtering" editable="true">
        <module label="SFilter v1" id="P_Filter" editableInJob="true">
            <description>This module filters reads based on a minimum subread length, polymerase read quality and polymerase read length.</description>
            <param name="minSubReadLength" label="Minimum Subread Length">
                <value>50</value>
                <title>Subreads shorter than this value (in base pairs) are filtered out and excluded from analysis.</title>
                <input type="text" size="3"/>
                <rule type="number" min="0.0" message="Value must be a positive integer"/>
            </param>
            <param name="readScore" label="Minimum Polymerase Read Quality">
                <value>75</value>
                <title>Polymerase reads with lower quality than this value are filtered out and excluded from analysis.</title>
                <input type="text"/>
                <rule type="number" min="0.0" message="Value must be between 0 and 100" max="100.0"/>
            </param>
            <param name="minLength" label="Minimum Polymerase Read Length">
                <value>50</value>
                <title>Polymerase reads shorter than this value (in base pairs) are filtered out and excluded from analysis.</title>
                <input type="text" size="3"/>
                <rule type="number" min="0.0" message="Value must be a positive integer"/>
            </param>
        </module>
        <module label="SFilter Reports v1" id="P_FilterReports" editableInJob="false"/>
    </moduleStage>
    <moduleStage name="barcode" editable="true"/>
    <fileName>RS_Subreads.1.xml</fileName>
</smrtpipeSettings>' . "\n";
close OF;

my $cmd1 = "smrtpipe.py -D TMP=$path/mytemp -D SHARED_DIR=$path/mytemp --output=$path --params=$tsettings xml:$txml";
print "$cmd1\n";
open(INF,"$cmd1 |") or die; 
while(my $line = <INF>) {
  chomp($line);
  print "$line\n";
}
close INF;
`cp $path/data/filtered_subreads.fasta $basefile.fasta`;
`cp $path/data/filtered_subreads.fastq $basefile.fastq`;
`rm -r $path`;
print "finished everything.. deleting folder $path and exiting.\n";
