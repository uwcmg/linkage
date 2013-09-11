#!/usr/bin/env perl
#
# Description:
#
#
#
# Created by Jessica Chong on 2013-06-01.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($rsIDcm2bp, $merlintbl, $cutofflod, $outputfile, $usechrprefix, $allchr, $help);

GetOptions(
	'cm2bp=s' => \$rsIDcm2bp,
	'merlintbl=s' => \$merlintbl,
	'allchr=s' => \$allchr,
	'cutofflod=f' => \$cutofflod,
	'usechr=s' => \$usechrprefix,
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>2, -exitval=>1) if $help;
if (!defined $rsIDcm2bp) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --cm2bp not defined.\n");
} elsif (!defined $merlintbl) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --merlintbl not defined.\n");
} elsif (!defined $cutofflod) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --cutofflod not defined.\n");
} elsif (!defined $usechrprefix) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --usechr not defined.\n");
} elsif (!defined $outputfile) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} elsif (!defined $allchr) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --allchr not defined\n");
}

my $maxlodseen = 0;
my $maxlodchr = 0;

open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
print $output_handle "#chr\tstartbp\tendbp\tname\tmaxLOD\n";

if ($allchr == 1 || $allchr =~ /TRUE/i || $allchr eq 'T' || $allchr eq 'Y' || $allchr =~ /YES/i) {
	$merlintbl =~ /([\.\/\w]+)\.chr(\w+)-([\.\w]+)/;
	my $merlinprefix = "$1.chr";			# et.dominant.chr10-parametric.tbl
	my $merlinsuffix = "-$3";
	
	for (my $chr=1; $chr<=22; $chr++) {
		getlinkageregions($rsIDcm2bp, "$merlinprefix$chr$merlinsuffix", $output_handle);
	}
} else {
	getlinkageregions($rsIDcm2bp, $merlintbl, $output_handle);
}

close $output_handle;


print "\nMax LOD seen is $maxlodseen on chromosome $maxlodchr\n";






sub getlinkageregions {
	my ($rsIDcm2bp, $merlintbl, $output_handle) = @_;
	my %map;
	open (my $plink_handle, "$rsIDcm2bp") or die "Cannot read $rsIDcm2bp: $!.\n";
	while ( <$plink_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $snpname, $cm, $pos) = split("\t", $_);
		$map{$snpname} = $pos;
	}
	close $plink_handle;


	my ($start, $stop, $inregion, $regionmaxlod, $currchr, $prevpos) = ((0) x 6);

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
		if (!defined $map{$snpname}) {
			print "$snpname doesn't exist in $chr??\n";
			exit;
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
				$start = $prevpos;
				$stop = $pos;
				$regionmaxlod = $lod;
				$inregion = 1;
			}
		} elsif ($lod < $cutofflod) {
			if ($inregion == 1) {
				$stop = $pos;
				print $output_handle "$printchr\t$start\t$stop\t$line[$analysiscol] $chr:$start-$stop\t$regionmaxlod\n";
				($start, $stop, $inregion, $regionmaxlod) = (0, 0, 0, 0);
			}
		}
		$prevpos = $pos;
	}
	close $merlintbl_handle;
}








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

