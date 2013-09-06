#!/usr/bin/env perl
#
# Description:
#
# Usage: perl untitled.pl
#
#
# Created by Jessica on 2011-11-21


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

GetOptions(
	'in=s' => \$inputfile, 
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;


for (my $i=1; $i<=22; $i++) {
	
	open (OUT, ">merlinfiles/CLBA.final.30k.merlin.chr$i.freq");
	open (FILE, "freqfiles/CLBA.final.30k.premerlin.chr$i.frq") or die "Cannot open  file.\n";
	<FILE>;
	while ( <FILE> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my @line = split (/\s+/, $_);
		my $major = (1-$line[5]);
		my $minor = $line[5];
		print OUT "M $line[2]\nF $minor\nF $major\n";
	}
	close FILE;
	close OUT;
	
}




################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


make_merline_freqfiles.pl - 


=head1 SYNOPSIS


perl B<make_merline_freqfiles.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--in> I<input file>	

	input file

=item B<--out> I<output file>

	name of output file

=back


=head1 DOCUMENTATION


xx


=head1 DESCRIPTION


This script will 


=head1 CAVEATS





=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
