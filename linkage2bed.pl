#!/usr/bin/env perl
#
# Description:
#
#
#
# Created by Jessica Chong on 20xx-xx-xx.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($plinkmap, $merlintbl, $cutofflod, $outputfile, $usechrprefix, $help);

GetOptions(
	'plinkmap=s' => \$plinkmap,
	'merlintbl=s' => \$merlintbl,
	'cutofflod=f' => \$cutofflod,
	'chrprefix=s' => \$usechrprefix,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;
if (!defined $plinkmap) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --plinkmap not defined.\n");
} elsif (!defined $merlintbl) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --merlintbl not defined.\n");
} elsif (!defined $cutofflod) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --cutofflod not defined.\n");
} elsif (!defined $usechrprefix) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --usechrprefix not defined.\n");
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} 


my %map;
open (my $plink_handle, "$plinkmap") or die "Cannot read $plinkmap: $!.\n";
while ( <$plink_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($chr, $snpname, $cm, $pos) = split("\t", $_);
	$map{$snpname} = $pos;
}
close $plink_handle;



my $maxlodseen = 0;
my $maxlodchr = 0;
my ($start, $stop, $inregion, $regionmaxlod, $currchr) = ((0) x 5);

open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
print $output_handle "#chr\tstart\tend\tname\tLOD\n";
open (my $merlintbl_handle, "$merlintbl") or die "Cannot read $merlintbl: $!.\n";
my $headerline = <$merlintbl_handle>;
my @header = split("\t", $headerline);
my ($lodcolumn, $analysiscol);
for (my $i=0; $i<=$#header; $i++) {
	if ($header[$i] eq 'LOD') {
		$lodcolumn = $i;
	}
	if ($header[$i] eq 'MODEL' || $header[$i] eq 'ANALYSIS') {
		$analysiscol = $i;
	}
}
while ( <$merlintbl_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split("\t", $_);
	my ($chr, $gmap, $snpname) = @line[0..2];
	if ($chr eq 'na') {
		next;
	}
	my $pos = $map{$snpname};
	my $lod = $line[$lodcolumn];
	if ($lod > $maxlodseen) {
		$maxlodseen = $lod;
		$maxlodchr = $chr;
	}
	
	my $printchr = $chr;
	if ($usechrprefix eq '1' || $usechrprefix =~ /T/i || $usechrprefix =~ /Y/i) {
		$printchr = "chr$chr";
	}
	
	if ($chr ne $currchr) {
		if ($inregion == 1) {
			print $output_handle "$printchr\t$start\t$stop\t$chr:$start-$stop\t$regionmaxlod\n";
		}
		($start, $stop, $inregion, $regionmaxlod) = (0, 0, 0, 0);
		$currchr = $chr;
	} 
	
	if ($lod >= $cutofflod) {
		if ($inregion == 1) {
			$stop = $pos;
			if ($lod > $regionmaxlod) {
				$regionmaxlod = $lod;
			}
		} else {
			$start = $pos;
			$stop = $pos;
			$regionmaxlod = $lod;
			$inregion = 1;
		}
	} elsif ($lod < $cutofflod) {
		if ($inregion == 1) {
			print $output_handle "$printchr\t$start\t$stop\t$line[$analysiscol] $chr:$start-$stop\t$regionmaxlod\n";
			($start, $stop, $inregion, $regionmaxlod) = (0, 0, 0, 0);
		}
	}
}
close $merlintbl_handle;
close $output_handle;


print "\nMax LOD seen is $maxlodseen on chromosome $maxlodchr\n";



################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


xxx.pl - 


=head1 SYNOPSIS


perl B<xxxx.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--in> F<input file>	

	input file

=item B<--out> F<output file>

	name of output file

=item B<--help> I<help>

	print documentation

=back


=head1 DOCUMENTATION


xx


=head1 DESCRIPTION


xxxxxxxx


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut

