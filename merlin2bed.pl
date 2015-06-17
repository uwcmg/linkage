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

my ($rsIDcm2bp, $merlintbl, $cutofflod, $outprefix, $usechrprefix, $allchr, $help);

GetOptions(
	'cm2bp=s' => \$rsIDcm2bp,
	'merlintbl=s' => \$merlintbl,
	'allchr=s' => \$allchr,
	'cutofflod=f' => \$cutofflod,
	'usechr=s' => \$usechrprefix,
	'outprefix=s' => \$outprefix,
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
} elsif (!defined $outprefix) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --out not defined\n");
} elsif (!defined $allchr) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --allchr not defined\n");
}

my $maxlodseen = 0;
my $maxlodchr = 0;

open (my $output_handle, ">", "$outprefix.bed") or die "Cannot write to $outprefix.bed: $!.\n";
open (my $output_extended_handle, ">", "$outprefix.detail.tsv") or die "Cannot write to $outprefix.detail.tsv: $!.\n";
print $output_handle "#chr\tstartbp\tendbp\tname\tmaxLOD\n";
print $output_extended_handle "#chr\tstartbp\tendbp\tname\tmaxLOD\tstartLOD\tendLOD\tstartSNP\tendSNP\tnSNPs\n";

if ($allchr =~ /T/i || $allchr =~ /Y/i || $allchr =~ /1/) {
	$merlintbl =~ /([\.\/\w]+)\.chr(\w+)-([\.\w]+)/;
	my $merlinprefix = "$1.chr";			# pheno.dominant.chr10-parametric.tbl
	my $merlinsuffix = "-$3";
	
	foreach my $chr ((1..22, "X")) {
		getlinkageregions($rsIDcm2bp, "$merlinprefix$chr$merlinsuffix", $output_handle, $output_extended_handle);
	}
} else {
	getlinkageregions($rsIDcm2bp, $merlintbl, $output_handle, $output_extended_handle);
}

close $output_handle;
close $output_extended_handle;


print "\nMax LOD seen is $maxlodseen on chromosome $maxlodchr\n";
print "Created $outprefix.bed and $outprefix.detail.tsv to summarize results of linkage analysis.\n";




sub getlinkageregions {
	my ($rsIDcm2bp, $merlintbl, $output_handle, $output_extended_handle) = @_;
	my %map;
	open (my $plink_handle, "$rsIDcm2bp") or die "Cannot read $rsIDcm2bp: $!.\n";
	while ( <$plink_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $snpname, $cm, $pos) = split("\t", $_);
		$map{$snpname} = $pos;
	}
	close $plink_handle;


	my ($inregion, $regionmaxlod, $currchr, $nsnps) = ((0) x 4);
	my @bpinfo = ((0) x 3);			# start bp, end bp, prev bp
	my @snpinfo = ((0) x 3);		# start snp, end snp, prevsnp
	my @lodinfo = ((0) x 3);		# start lod, end lod, prevlod
	
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
		if ($line[3] =~ /\[Pairs\]/) {		# if this is a non-parametric file, it'll contain values for ALL and PAIRS statistics
			last;
		}
		if ($chr eq 'na') {
			next;
		}
		if (!defined $map{$snpname}) {
			die "$snpname doesn't exist in $chr??\n";
		}
		my $pos = $map{$snpname};
		my $lod = $line[$lodcolumn];
		if ($lod > $maxlodseen) {
			$maxlodseen = $lod;
			$maxlodchr = $chr;
		}
	
		my $printchr = $chr;
		if ($usechrprefix =~ /1/ || $usechrprefix =~ /T/i || $usechrprefix =~ /Y/i) {
			$printchr = "chr$chr";
		}
	
		if ($chr ne $currchr) {
			if ($inregion == 1 && $nsnps > 1) {
				print "ending peak,newchr: $snpname / $chr:$pos ($nsnps) = $lod\n";
				print $output_handle "$printchr\t$bpinfo[0]\t$bpinfo[1]\t$chr:$bpinfo[0]-$bpinfo[1]\t$regionmaxlod\n";
				print $output_extended_handle "$printchr\t$bpinfo[0]\t$bpinfo[1]\t$chr:$bpinfo[0]-$bpinfo[1]\t$regionmaxlod\t$lodinfo[0]\t$lodinfo[1]\t$snpinfo[0]\t$snpinfo[1]\t$nsnps\n";
			}
			($inregion, $regionmaxlod, $nsnps) = ((0) x 3);
			@bpinfo = ((0) x 3);			# start bp, end bp, prev bp
			@snpinfo = ((0) x 3);		# start snp, end snp, prevsnp
			@lodinfo = ((0) x 3);		# start lod, end lod, prevlod
			$currchr = $chr;
		} 
		
		if ($lod >= $cutofflod) {
			if ($inregion == 1) {
				# print "already in peak: $snpname / $chr:$pos ($nsnps) = $lod\n";
				# already in a linkage peak, update the ends to match the current locus
				$nsnps++;
				$bpinfo[1] = $pos;
				$lodinfo[1] = $lod;
				$snpinfo[1] = $snpname;
				if ($lod > $regionmaxlod) {
					$regionmaxlod = $lod;
				}
			} else {
				# print "\nstarting peak: $snpname / $chr:$pos ($nsnps) = $lod\n";
				# starting a NEW linkage peak
				$inregion = 1;
				$regionmaxlod = $lod;
				$nsnps = 1;
				# the start should be the info for the previous locus
				$bpinfo[0] = $bpinfo[2];
				$lodinfo[0] = $lodinfo[2];
				$snpinfo[0] = $snpinfo[2];
				# the ends should be the info for the current locus
				$bpinfo[1] = $pos;
				$lodinfo[1] = $lod;
				$snpinfo[1] = $snpname;
			}
		} elsif ($lod < $cutofflod) {
			if ($inregion == 1) {
				# print "ending peak: $snpname / $chr:$pos ($nsnps) = $lod\n";
				# have left a linkage peak
				$bpinfo[1] = $pos;
				$snpinfo[1] = $snpname;
				$lodinfo[1] = $lod;
				if ($nsnps > 1) {
					$nsnps++;
					print $output_handle "$printchr\t$bpinfo[0]\t$bpinfo[1]\t$chr:$bpinfo[0]-$bpinfo[1]\t$regionmaxlod\n";
					print $output_extended_handle "$printchr\t$bpinfo[0]\t$bpinfo[1]\t$chr:$bpinfo[0]-$bpinfo[1]\t$regionmaxlod\t$lodinfo[0]\t$lodinfo[1]\t$snpinfo[0]\t$snpinfo[1]\t$nsnps\n";
				}
				($inregion, $regionmaxlod, $nsnps) = ((0) x 3);
				@bpinfo = ((0) x 3);			# start bp, end bp, prev bp
				@snpinfo = ((0) x 3);		# start snp, end snp, prevsnp
				@lodinfo = ((0) x 3);		# start lod, end lod, prevlod
			}
		}
		
		if ($inregion == 1 && $nsnps > 1 && eof) {
			# print "ending peak, EOF\n";
			print $output_handle "$printchr\t$bpinfo[0]\t$bpinfo[1]\t$chr:$bpinfo[0]-$bpinfo[1]\t$regionmaxlod\n";
			print $output_extended_handle "$printchr\t$bpinfo[0]\t$bpinfo[1]\t$chr:$bpinfo[0]-$bpinfo[1]\t$regionmaxlod\t$lodinfo[0]\t$lodinfo[1]\t$snpinfo[0]\t$snpinfo[1]\t$nsnps\n";
		}
		
		# save this info in case we start a new peak here
		$bpinfo[2] = $pos;
		$snpinfo[2] = $snpname;
		$lodinfo[2] = $lod;
	}
	close $merlintbl_handle;


}








################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


merlin2bed.pl - Given results from Merlin (using --markerNames --tabulate), output a BED file with coordinates of regions meeting a LOD score cutoff and a BED-like file (with extra columns) with more information on the regions.


=head1 SYNOPSIS


perl B<merlin2bed.pl> I<[options]>


=head1 ARGUMENTS


=over 4


=item B<--cm2bp> F<*merlin.cm2bp.map>	

	PLINK map format to provide the physical positions of variants used in linkage analysis

=item B<--merlintbl> I<.tbl merlin file>

	path to merlin results file (produced when you use --tabulate --markerNames)
	
=item B<--allchr> I<TRUE|FALSE>

	should this script look across all chromosomes (1:22) with same filename pattern as the --merlintbl argument or just use the one input file

=item B<--cutofflod> I<number>

	minimum LOD score for a marker to be considered part of a linkage peak (peaks will include all markers with LOD>=cutoff plus flanking markers)

=item B<--usechr> I<TRUE|FALSE>

	use the "chr22" notation instead of "22" in the output BED file

=item B<--outprefix> F<filename>

	prefix for files produced by this script

=item B<--help> I<help>

	print documentation

=back


=head1 NOTES


If plink2merlin.pl was used to generate the input files for Merlin, then there should be a file named <phenotype>.merlin.cm2bp.map in the same directory as the Merlin files and this file can be used for the --cm2bp argument.

If you have the nextgen data in VCF or BED format, you can Bedtools to get the variants within the regions in the output from this script using a command like:

intersectBed -a myphenotype.nextgenvariants.tsv -b linkageresults.bed -header > myphenotype.overlaplinkage.tsv

intersectBed -a myphenotype.exome.vcf -b linkageresults.bed -header > myphenotype.overlaplinkage.tsv


=head1 EXAMPLES


	perl merlin2bed.pl
		--cm2bp myphenotype.merlin.map 
		--merlintbl myphenotype.dominant.chr22-parametric.tbl 
		--cutofflod 1 
		--usechr T 
		--outprefix temp 
		--allchr T


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)

