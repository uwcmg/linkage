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


my ($freqdir, $mapdir, $outdir, $minfreq, $maxfreq, $chip, $mincMdist, $help);

GetOptions(
	'freqdir=s' => \$freqdir, 
	'mapdir=s' => \$mapdir, 
	'outdir=s' => \$outdir,
	'minfreq=f' => \$minfreq,
	'maxfreq=f' => \$maxfreq,
	'cMdist=f' => \$mincMdist,
	'chip=s' => \$chip,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>2) if $help;

if (!defined $freqdir) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --freqdir not defined.\n")
} elsif (!defined $mapdir) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --mapdir not defined\n");
} elsif (!defined $outdir) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --outdir not defined\n");
} elsif (!defined $minfreq) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --minfreq not defined\n");
} elsif (!defined $maxfreq) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --maxfreq not defined\n");
} elsif (!defined $chip) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --chip not defined\n");
} elsif (!defined $mincMdist) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --cMdist not defined\n");
}    
  
$freqdir = remove_trailing_slash_dir($freqdir);
$mapdir = remove_trailing_slash_dir($mapdir);
$outdir = remove_trailing_slash_dir($outdir);

my $logfile = "$chip.highhet.log";

open (my $log_handle, ">", "$outdir/$logfile") or die "Cannot write to $outdir/$logfile: $!.\n";
my $totalpassfreqvar = 0;
my $totalgridvar = 0;
my $totalvar = 0;
foreach my $chr ((1..22)) {
	print $log_handle "Alt allele frequency in each population must be between:  $minfreq to $maxfreq\n";
	print $log_handle "Grid markers must be spaced at least $mincMdist cM apart\n\n";

	print $log_handle "Selecting high heterozygosity SNPs on chr $chr for $chip\n";
	print "Selecting high heterozygosity SNPs on chr $chr for $chip\n";
	# read in sex-averaged cM positions
	my %geneticmap;
	open (my $map_handle, "$mapdir/chr$chr.$chip.map") or die "Cannot read $mapdir/chr$chr.$chip.map: $!.\n";
	while ( <$map_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($SNPid, $rsID, $bp37, $sex_avg_kos, $female_kos, $male_kos, $sex_avg_hald, $female_hald, $male_hald) = split("\t", $_);
		$geneticmap{$rsID} = $sex_avg_hald;
		$geneticmap{$SNPid} = $sex_avg_hald;
	}
	close $map_handle;
	
	# read through to get frequencies
	my $outputfile = "chr$chr.$chip.highhet.txt";
	open (my $output_handle, ">", "$outdir/$outputfile") or die "Cannot write to $outdir/$outputfile: $!.\n";
	print $output_handle "SNPid\trsID\tbp37\tref\talt\tAFR\tAMR\tASN\tEUR\tOVERALL\tsex_avg_hald\tpasscMdist\n";
	
	open (my $gridallchr_handle, ">", "$outdir/chr$chr.$chip.grid.snplist") or die "Cannot write to chr$chr.$chip.grid.snplist: $!.\n";

	my $prevcM = -100;
	my $countpassfreqvar = 0;
	my $countgridvar = 0;
	open (my $freq_handle, "$freqdir/chr$chr.$chip.freq") or die "Cannot read $freqdir/chr$chr.$chip.freq: $!.\n";
	while ( <$freq_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my @line = split("\t", $_);
		my ($SNPid, $rsID, $b37, $ref, $alt, @popfreqs) = split("\t", $_);  	# order of population allele freqs: $AFR, $AMR, $ASN, $EUR, $OVERALL 
		if (ishighhet(\@popfreqs, $minfreq, $maxfreq)) {
			if (defined $geneticmap{$rsID}) {
				print $output_handle join("\t", @line)."\t$geneticmap{$rsID}\t";
				if ($geneticmap{$rsID} >= ($prevcM + $mincMdist)) {
					print $output_handle "Y\n";
					$prevcM = $geneticmap{$rsID};
					$countgridvar++;
					print $gridallchr_handle "$rsID\n";
				} else {
					print $output_handle "N\n";
				}
				$countpassfreqvar++;
			} else {
				print $log_handle "\t$SNPid (rsID=$rsID) doesn't have a genetic map position on chr$chr\n";
			}
		} 
		$totalvar++;
	}
	close $freq_handle;
	close $output_handle;
	close $gridallchr_handle;
	
	print $log_handle "$countpassfreqvar markers on chr$chr had all population frequencies >=$minfreq and <=$maxfreq and $countgridvar pass genetic intermarker distance >= $mincMdist cM\n";
	$totalpassfreqvar += $countpassfreqvar;
	$totalgridvar += $countgridvar;
}

print $log_handle "\n";
print $log_handle "$totalpassfreqvar markers out of $totalvar had all population frequencies >=$minfreq and <=$maxfreq, resulting in total $totalgridvar variants for the grid\n\n";
close $log_handle;

print "Selected $totalgridvar variants for the grid\n\n";
print "View $logfile for details\n";






sub ishighhet {
	my ($popfreq_ref, $minfreq, $maxfreq) = @_;
	my $ishighhet = 1;
	my @subpopfreqs = @{$popfreq_ref};
	pop(@subpopfreqs);
	
	foreach my $freq (@subpopfreqs) {
		if ($freq < $minfreq || $freq > $maxfreq) {
			$ishighhet = 0;
		}
	}
	
	return $ishighhet;
}

sub remove_trailing_slash_dir {
	my $dirname = $_[0];
	$dirname =~ s/[\/|\\]$//;
	if (! -e $dirname) {
		die "$dirname does not exist\n";
	}
	return $dirname;
}



################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


get_highhet_variants.pl - Given freq and map files (Liz's format) for a given chip, generate a list of markers to be used for a linkage grid panel.


=head1 SYNOPSIS


perl B<get_highhet_variants.pl> I<[options]>


=head1 ARGUMENTS

=over 4

=item B<--freqdir> F<input file>	

	path to directory containing chr$chr.$chip.freq files with population allele frequencies

=item B<--mapdir> F<input file>	

	path to directory containing chr$chr.$chip.map files with haldane distsances

=item B<--outdir> F<input file>	

	directory to save files produced by this script

=item B<--minfreq> F<input file>	

	required minimum allele frequency in any population (if lower, do not use marker in grid)

=item B<--maxfreq> F<input file>	

	required maximum allele frequency in any population (if higher, do not use marker in grid)

=item B<--cMdist> F<input file>	

	minimum haldane distance (cM) between markers in grid

=item B<--chip> F<ExomeChip|CytoChip>	

	name of genotyping chip

=item B<--help> I<help>

	print documentation

=back


=head1 DESCRIPTION


For each variant in the freq files, check if the alt allele frequency is between (minfreq) and (maxfreq) for all populations in the file, if the variant passes this check, make sure it has a sex-averaged Haldane map position available, then check if marker is at least (cMdist) centimorgans away from the previous marker.  If so, add marker to the linkage grid list.  If not, continue onto next marker.  


=head1 EXAMPLES


	perl get_highhet_variants.pl 
		--freqdir ../../CytoChipComplete/freqs/ 
		--mapdir ../../CytoChipComplete/maps/ 
		--outdir . 
		--minfreq 0.3 
		--maxfreq 0.7 
		--cMdist 0.5 
		--chip CytoChip


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
