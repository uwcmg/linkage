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


my ($freqdir, $mapdir, $outdir, $minfreq, $maxfreq, $chip, $help);

GetOptions(
	'freqdir=s' => \$freqdir, 
	'mapdir=s' => \$mapdir, 
	'outdir=s' => \$outdir,
	'minfreq=f' => \$minfreq,
	'maxfreq=f' => \$maxfreq,
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
				if ($geneticmap{$rsID} >= ($prevcM + 0.7)) {
					print $output_handle "Y\n";
					$prevcM = $geneticmap{$rsID};
					$countgridvar++;
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
	
	print $log_handle "$countpassfreqvar markers on chr$chr had all population frequencies >=$minfreq and <=$maxfreq and $countgridvar pass genetic intermarker distance >= 0.7 cM\n";
	$totalpassfreqvar += $countpassfreqvar;
	$totalgridvar += $countgridvar;
}


print $log_handle "\n";
print $log_handle "Total $totalpassfreqvar markers out of $totalvar had all population frequencies >=$minfreq and <=$maxfreq and total $totalgridvar variants in the grid\n\n";

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


=head1 FILES


xx


=head1 EXAMPLES


xxxxxxxx


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
