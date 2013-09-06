#!perl
#
# Description:
#
# Usage: perl add_cM_pos.pl <Merlin format MAPFILE>
#
#
# Created by Jessica on 2010-03-22

use strict;
use warnings;

if (@ARGV != 1) {
	print "Usage: perl add_cM_pos.pl <Merlin MAPFILE without header> > output.merlin.map\n";
	exit;
}

my $mapfile = $ARGV[0];


my %map;
open (FILE, "affygeneticmap/500k.geneticmap") or die "Cannot open  file.\n";
while ( <FILE> ){
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split ("\t", $_);
	$map{$line[0]} = $line[1];
}
close FILE;


open (FILE, "$mapfile") or die "Cannot open $mapfile file.\n";
print "CHROMOSOME MARKER LOCATION\n";
while ( <FILE> ){
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split (" ", $_);
	my $markername = $line[1];
	if (!exists $map{$markername}) {
		print "$markername doens't exist in map\n";
		die;
	} else {
		print "$line[0] $line[1] $map{$markername}\n";
	}
}
close FILE;


