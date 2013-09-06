#!perl
#
# Description:
#
# Usage: perl getLod.pl
#
#
# Created by Jessica on 2010-05-17

use strict;
use warnings;

if (@ARGV != 4) {
	print "\nUsage: perl getLod.pl PHENOTYPE_NAME FOLDER CUTOFF_FOR_LODs chr/combined\n\n";
	exit;
}

my $pheno = $ARGV[0];
my $targetfolder = $ARGV[1];
$targetfolder =~ s/\///;
my $cutofflod = $ARGV[2];
my $files = $ARGV[3];

my $maxlodseen = 0;
my $maxlodchr = 0;

if ($files eq "chr") {
	print "Getting LOD scores >= $cutofflod from all chromosomes\nSaving to file: $targetfolder/$pheno.lodsummary.txt\n";
	my $outputfile = "$targetfolder/$pheno.lodsummary.txt";
	my $model = `cat $targetfolder/$pheno.model`;
	my $header = `head -1 "$targetfolder/$pheno.chr1.pl-parametric.tbl"`;
	open (OUTPUT, ">$outputfile") or die "Cannot write to $outputfile .\n";
	print OUTPUT "$model\n\n";
	print OUTPUT "$header";
	for (my $i=1; $i<=22; $i++) {
		my $resulttable = "$targetfolder/$pheno.chr$i.pl-parametric.tbl";
		open (FILE, "$resulttable") or die "Cannot open $resulttable file.\n";
		<FILE>;
		while ( <FILE> ){
			$_ =~ s/\s+$//;					# Remove line endings
			my @line = split ("\t", $_);
			my $lod = $line[4];
			if ($lod > $maxlodseen) {
				$maxlodseen = $lod;
				$maxlodchr = $line[0];
			}
			if ($lod >= $cutofflod) {
				print OUTPUT join("\t", @line)."\n";
			}
		}
		close FILE;
	}
	close OUTPUT;
	print "\nMax LOD seen is $maxlodseen on chromosome $maxlodchr\n";
}

if ($files eq "combined") {
	print "Getting LOD scores >= $cutofflod from all chromosomes\nSaving to file: $targetfolder/$pheno.lodsummary.txt\n";
	my $outputfile = "$targetfolder/$pheno.lodsummary.txt";
	my $header = `head -1 "$targetfolder/$pheno-parametric.tbl"`;
	my $model = `cat $targetfolder/$pheno.model`;
	open (OUTPUT, ">$outputfile") or die "Cannot write to $outputfile .\n";
	print OUTPUT "$model\n\n";
	print OUTPUT "$header";
	my $resulttable = "$targetfolder/$pheno-parametric.tbl";
	open (FILE, "$resulttable") or die "Cannot open $resulttable file.\n";
	<FILE>;
	while ( <FILE> ){
		$_ =~ s/\s+$//;					# Remove line endings
		my @line = split ("\t", $_);
		my $lod = $line[4];
		if ($lod > $maxlodseen) {
			$maxlodseen = $lod;
			$maxlodchr = $line[0];
		}
		if ($lod >= $cutofflod) {
			print OUTPUT join("\t", @line)."\n";
		}
	}
	close FILE;
	close OUTPUT;
	print "\nMax LOD seen is $maxlodseen on chromosome $maxlodchr\n";
}

