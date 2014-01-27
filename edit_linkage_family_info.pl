#!perl
#
# Description:
#
# Usage: perl edit_linkage_family_info.pl
#
#
# Created by Jessica on 2010-05-11

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my $editsfile;
my $pedfile;
my $outputfile;
my $help;

GetOptions(
	'edits=s' => \$editsfile, 
	'pedfile=s' => \$pedfile, 
	'out=s' => \$outputfile,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>1, -exitval=>1) if $help;



my %info;
my %fakesubjects;
my %removesubjects;
open (FILE, "$editsfile") or die "Cannot open $editsfile file.\n";
while ( <FILE> ){
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split (/\s+/, $_);
	my $subjectid = $line[1];
	if ($subjectid =~ s/#//) {
		$line[1] =~ s/#//;
		$fakesubjects{$subjectid} = 1;
	}
	if ($subjectid =~ s/\!//) {
		$line[1] =~ s/\!//;
		$removesubjects{$subjectid} = 1;
	}
	@{$info{$subjectid}} = @line;
}
close FILE;


open (my $output_handle, ">", "$outputfile") or die "Cannot write to $outputfile: $!.\n";
open (FILE, "$pedfile") or die "Cannot open $pedfile file.\n";
my $countalleles;
while ( <FILE> ){
	$_ =~ s/\s+$//;					# Remove line endings
	my @line = split (" ", $_);
	my $subjectid = $line[1];
	if (defined @{$info{$subjectid}} && !defined $removesubjects{$subjectid}) {
		print $output_handle join(' ', @{$info{$subjectid}})." @line[6..$#line]\n";	
	} 
	# elsif (!defined $removesubjects{$subjectid}) {
		# print $output_handle join(' ', @line)."\n";	
	# }
	$countalleles = $#line + 1 - 6;			
}
close FILE;


# $emptygenotypes =~ s/A|T|G|C|1|2/0/g;
my $emptygenotypes = join(' ', ((0) x $countalleles));
# Add in fake genotypes for people who need to be added artificially
foreach my $subjectid (keys %fakesubjects) {
	print $output_handle join(' ', @{$info{$subjectid}})." $emptygenotypes\n";
}
close $output_handle;



################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


edit_linkage_family_info.pl - 


=head1 SYNOPSIS


perl B<edit_linkage_family_info.pl> I<[options]>


=head1 ARGUMENTS


=over 4

=item B<--edits> I<input file>	

	contains familyid subject father mother sex affection_status  :  cut -f1-6 -d' ' FILE > output.txt  NOTE: add a # sign before the subject id number if the person is being added with empty genotypes


=item B<--pedfile> I<input file>	

	PLINK .ped file containing genotypes and family information


=item B<--out> I<output file>

	name of output file

=back


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut

