#!/usr/bin/env perl
#
# Description: Given the PLINK-formatted out by Bead/GenomeStudio, rename variants to rsIDs, get genetic map coordinates (haldane) from Liz's files, update maps, update family and parent information, zero out Mendelian errors, create per-chromosome PLINK files, convert to Merlin format, set up script to submit linkage jobs.
#
#
#
# Created by Jessica Chong on 2013-09-04.

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Copy;

my $mendeliandir = '/net/grc/vol1/mendelian_projects/mendelian_analysis/references';
my ($rawgenodir, $genotypechip, $refpop, $pheno, $model, $familyedits, $interimdir, $outdir, $help);

GetOptions(
	'rawgenodir=s' => \$rawgenodir, 
	'chip=s' => \$genotypechip, 
	'refpop=s' => \$refpop,
	'pheno=s' => \$pheno,
	'model=s' => \$model, 
	'familyedits:s' => \$familyedits, 
	'interimdir=s' => \$interimdir,
	'outdir=s' => \$outdir,
	'help|?' => \$help,
) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose=>2, -exitval=>1) if $help;

if (!defined $rawgenodir) {
	 pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --rawgenodir not defined.\n")
} elsif (!defined $genotypechip) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --chip not defined\n");
} elsif (!defined $refpop) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --refpop not defined\n");
} elsif (!defined $pheno) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --pheno not defined\n");
} elsif (!defined $interimdir) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --interimdir not defined\n");
} elsif (!defined $outdir) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --outdir not defined\n");
} elsif (!defined $model) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --model not defined\n");
} 


$rawgenodir = remove_trailing_slash_dir($rawgenodir);
$outdir = remove_trailing_slash_dir($outdir);
$interimdir = remove_trailing_slash_dir($interimdir);

my ($updatecMfile, $refdatadir);


print "\n";
print "... using $genotypechip linkage markers and frequencies\n";	
if ($genotypechip =~ /ExomeChip/i) {
	$refdatadir = "$mendeliandir/ExomeChipLinkage";
	if (system("cut -f1,7 $mendeliandir/ExomeChipCompleteMaps/chr*.ExomeChip.map | sort -k1 | uniq > $interimdir/allchr.ExomeChip.updatecM.txt") != 0) {
		die "Can't extract haldane values from $mendeliandir/ExomeChipCompleteMaps/chr*.ExomeChip.map: $?";
	}
	$updatecMfile = 'allchr.ExomeChip.updatecM.txt';
	if (system("cut -f1 $interimdir/allchr.ExomeChip.updatecM.txt > $interimdir/rsIDvariantsonly.snplist")) {
		die "Can't extract rsIDs with cM values from allchr.ExomeChip.updatecM.txt: $?";
	}  
} elsif ($genotypechip =~ /CytoChip/i) {
	die "Don't have genetic maps in $mendeliandir for $genotypechip\n";
}



print "... generating list of variants with rsIDs in GenomeStudio output and cleaning up names\n";
# read contents of the sample_qc/PLINK* directory and get names of map and ped files
opendir my($dir_handle), $rawgenodir or die "Couldn't open dir $mendeliandir: $!";
my @rawgenofiles = readdir $dir_handle;
closedir $dir_handle;
my ($rawgenomap, $rawgenoped, $rawgenostem);
foreach my $file (@rawgenofiles) {
	if ($file =~ /.ped$/) {
		$rawgenoped = $file;
	}
	if ($file =~ /.map$/) {
		$rawgenomap = $file;
	}
}
if (!defined $rawgenoped || !defined $rawgenomap) {
	die "Cannot find GenomeStudio PLINK-formatted ped/map files in $rawgenodir\n";
}
$rawgenostem = $rawgenomap;
$rawgenostem =~ s/.map//;



# generate list of variants with rsIDs that are in the raw genotype data and cleanup the name to be only the rsID itself
open (my $update_rsid_handle, ">", "$interimdir/update_rsIDs.txt") or die "Cannot write to $interimdir/update_rsIDs.txt: $!.\n";
open (my $raw_map_handle, "$rawgenodir/$rawgenomap") or die "Cannot read $rawgenodir/$rawgenomap: $!.\n";
while ( <$raw_map_handle> ) {
	$_ =~ s/\s+$//;					# Remove line endings
	my ($chr, $varname, $cM, $bp) = split("\t", $_);
	if ($varname =~ 'rs') {
		my $rsid = $varname;
		$rsid =~ s/exm-//;
		print $update_rsid_handle "$varname\t$rsid\n";
	}
}
close $raw_map_handle;
close $update_rsid_handle;


print "... updating PLINK files, extracting polymorphic (MAF>=0.01) sites with rsIDs and call rate >= 95%, updating genetic maps\n";
# update raw genotype files with rsIDs and create new PLINK files with only those variants
`plink --file $rawgenodir/$rawgenostem --update-map $interimdir/update_rsIDs.txt --update-name --make-bed --out $interimdir/$pheno.allvar`;
`plink --bfile $interimdir/$pheno.allvar --extract $interimdir/rsIDvariantsonly.snplist --make-bed --out $interimdir/$pheno.allrsIDs`;

# update PLINK files with genetic map information
`plink --bfile $interimdir/$pheno.allrsIDs --update-map $interimdir/$updatecMfile --update-cm --make-bed --out $interimdir/$pheno.haldane`;
# verify all variants have genetic map updated: "0 in data but not in [ /nfs/home/jxchong/lib/allchr.ExomeChip.updatecM.txt ]"

# do some basic QC
# minimum 95% genotyping rate --geno 0.05
# minimum MAF>= 0.01 in the entire set of subjects; mostly to exclude monomorphic markers --maf 0.01 --nonfounders
`plink --bfile $interimdir/$pheno.haldane --geno 0.05 --make-bed --out $interimdir/$pheno.callrate95`;
`plink --bfile $interimdir/$pheno.callrate95 --maf 0.01 --nonfounders --make-bed --out $interimdir/$pheno.polymorphic`;

# update family ID, phenotype, parental information based on files
`plink --bfile $interimdir/$pheno.polymorphic --update-ids $pheno.updateFID.txt --make-bed --out $interimdir/$pheno.updateFID`;
`plink --bfile $interimdir/$pheno.updateFID --update-parents $pheno.updateparents.txt --make-bed --out $interimdir/$pheno.updateparents`;

# flip alleles in genotype file to match the 1000 Genomes reference file
makefliplist("$pheno.updateparents.bim", $refdatadir);
`plink --bfile $interimdir/$pheno.updateparents --flip $interimdir/rsIDs.toflip.txt --exclude $interimdir/rsIDs.toexclude.txt --make-bed --out $pheno.flipped`;

# do mendelian error check and zero out all probematic genotypes
`plink --bfile $interimdir/$pheno.flipped --me 1 1 --set-me-missing --missing-phenotype 0 --recode --out $interimdir/$pheno.updateparents.me1-1`;
copy("$interimdir/$pheno.updateparents.me1-1.map", "$interimdir/$pheno.merlin.map") or die "Failed to copy $interimdir/$pheno.updateparents.me1-1.map to $interimdir/$pheno.merlin.map\n";


# add in any extra ancestors and exclude people if desired
if (defined $familyedits) {
	print "\tPedigree edits provided ($familyedits); making changes.\n";
	if (system("perl /nfs/home/jxchong/bin/linkage_scripts/edit_linkage_family_info.pl --edits $familyedits --pedfile $interimdir/$pheno.updateparents.me1-1.ped --out $interimdir/$pheno.merlin.ped") != 0) {
		die "Could not run: perl /nfs/home/jxchong/bin/linkage_scripts/edit_linkage_family_info.pl --edits $familyedits --pedfile $interimdir/$pheno.updateparents.me1-1.ped --out $interimdir/$pheno.merlin.ped: $?";
	}
} else {
	copy("$interimdir/$pheno.updateparents.me1-1.ped", "$interimdir/$pheno.merlin.ped") or die "Failed to copy $interimdir/$pheno.updateparents.me1-1.ped to $interimdir/$pheno.merlin.ped\n";
}


# create a default model file if one doesn't already exist
if (! -e "$outdir/$pheno.model") {
	if ($model =~ /dominant/i) {
		`echo '$pheno 0.001 0,0.95,1.0 Dominant_(q=0.001,_f=0,0.95,1.0)' | cat > $outdir/$pheno.model`;
	} elsif ($model =~ /recessive/i) {
		`echo '$pheno 0.005 0,0,1.0 Recessive_(q=0.005,_f=0,0,1.0)' | cat > $outdir/$pheno.model`;
	} else {
		`echo '$pheno 0.005 0,0,1.0 Generic_(q=0.005,_f=0,0,1.0)' | cat > $outdir/$pheno.model`;
	}
}


# Output linkage file per chromosome
print "... making linkage files\n";
for (my $i=1; $i<=22; $i++) {
	print "\tfor chromosome $i\n";
	# create ped file
	`plink --file $interimdir/$pheno.merlin --chr $i --missing-phenotype 0 --recode --out $interimdir/$pheno.merlin.chr$i`;
	move("$interimdir/$pheno.merlin.chr$i.ped", "$outdir/$pheno.chr$i.ped") or die "Failed to move $interimdir/$pheno.merlin.chr$i.ped to $outdir/$pheno.chr$i.ped\n";
	
	# create map file 
	`cut -f1-3 $interimdir/$pheno.merlin.chr$i.map > $outdir/$pheno.chr$i.map`;
	
	# create frequency file
	create_frqfile($refdatadir, $refpop, $i);
	# allele flipping to account for different strands

	# create .dat file
	`echo 'A $pheno' > $outdir/$pheno.chr$i.dat`;
	`perl -ane \'print \"M \$F[1]\\n\";\' $interimdir/$pheno.merlin.chr$i.map >> $outdir/$pheno.chr$i.dat`;
}

print "Merlin format files $pheno.chr*.map $pheno.chr*.ped $pheno.chr*.freq $pheno.chr*.dat $pheno.model are ready in $outdir.\n";
print "$pheno.model may need to be edited depending on your model of inheritance\n\n";



# create job submission script to run linkage
open (my $submit_handle, ">", "$outdir/sublinkage.sh") or die "Cannot write to $outdir/sublinkage.sh: $?.\n";
print $submit_handle "#$ -S /bin/bash\n";
print $submit_handle "#$ -t 1-22\n";
print $submit_handle "#$ -o .\n";
print $submit_handle "#$ -e .\n";
print $submit_handle "#$ -cwd\n\n";
print $submit_handle "merlin -d $pheno.chr\${SGE_TASK_ID}.dat -p $pheno.chr\${SGE_TASK_ID}.ped -m $pheno.chr\${SGE_TASK_ID}.map -f $pheno.chr\${SGE_TASK_ID}.freq --model $pheno.model --pdf > $pheno.chr\${SGE_TASK_ID}.merlin.out\n";

print $submit_handle "### in case you need to manually re-run linkage on one of the chromosomes:\n";
for (my $i=1; $i<=22; $i++) {
	print $submit_handle "## merlin -d $pheno.chr$i.dat -p $pheno.chr$i.ped -m $pheno.chr$i.map -f $pheno.chr$i.freq --model $pheno.model --pdf > $pheno.chr$i.merlin.out\n";
}
close $submit_handle;

print "To run linkage: qsub $outdir/sublinkage.sh\n";








################################################################################################################
############################################# Subfunctions #####################################################
################################################################################################################

sub remove_trailing_slash_dir {
	my $dirname = $_[0];
	$dirname =~ s/[\/|\\]$//;
	if (! -e $dirname) {
		die "$dirname does not exist\n";
	}
	return $dirname;
}

sub makefliplist {
	my ($bimfile, $refdatadir) = @_;	

	my %ref_alleles;
	for (my $chr=1; $chr<=22; $chr++) {
		open (my $ref_allele_handle, "$refdatadir/freqs/grid$chr.freqs") or die "Cannot read $refdatadir/freqs/grid$chr.freqs: $?.\n";
		while (<$ref_allele_handle>) {
			$_ =~ s/\s+$//;					# Remove line endings
			my ($SNPid, $rsID, $b37, $ref, $alt, @popfreqs) = split("\t", $_); 
			$ref_alleles{$rsID}{'ref'} = $ref;
			$ref_alleles{$rsID}{'alt'} = $alt;
		}
		close $ref_allele_handle;
	}
	
	open (my $fliplist_handle, ">", "$interimdir/rsIDs.toflip.txt") or die "Cannot write to $interimdir/rsIDs.toflip.txt: $?.\n";
	open (my $ambiguouslist_handle, ">", "$interimdir/rsIDs.toexclude.txt") or die "Cannot write to $interimdir/rsIDs.toexclude.txt: $?.\n";
	open (my $bim_handle, "$interimdir/$bimfile") or die "Cannot read $interimdir/$bimfile: $?.\n";
	while (<$bim_handle>) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $rsID, $cM, $b37, $a1, $a2) = split("\t", $_);
		my $action = needsflip($a1, $a2, $ref_alleles{$rsID}{'ref'}, $ref_alleles{$rsID}{'alt'});
		if ($action eq 'flip') {
			print $fliplist_handle "$rsID\n";
		} elsif ($action eq 'ambiguous' || $action eq 'monomorphic') {
			print $ambiguouslist_handle "$rsID\n";
		} elsif ($action eq 'weird') {
			print STDERR "$rsID on chr$chr is not ambiguous (AT or GC SNP) but doesn't need to be flipped either.\n";
			print $ambiguouslist_handle "$rsID\n";
		}
	}
	close $bim_handle;	
	close $ambiguouslist_handle;	
	close $fliplist_handle;
}

sub needsflip {
	my ($a1, $a2, $ref, $alt) = @_;
	my %flipallele = ('A'=>'T', 'C'=>'G', 'G'=>'C', 'T'=>'A', '0' => '0'); 			# switch the strand of the genotype

	my $action = 'weird';
	if ("$a1$a2" =~ '0') {
		$action = 'monomorphic';
	} elsif ("$a1$a2" eq 'AT' || "$a1$a2" eq 'TA' || "$a1$a2" eq 'GC' || "$a1$a2" eq 'CG') {
		$action = 'ambiguous';
	} elsif ("$a1$a2" =~ /$ref/ && "$a1$a2" =~ /$alt/) {
		$action = 'nochange';
	} elsif ("$a1$a2" =~ /$flipallele{$ref}/ && "$a1$a2" =~ /$flipallele{$alt}/) {
		$action = 'flip';
	}
	return $action;
}

sub create_frqfile {
	my $chr = $_[0];
	
	# determine column number for getting ref population's frequencies
	my $header = `head -1 $refdatadir/freqs/header.grid.freqs`;
	chomp $header;
	my @cols = split("\t", $header);
	my $refpopcol;
	for (my $colnum=5; $colnum<=$#cols; $colnum++) {
		if ($cols[$colnum] =~ /$refpop/i) {
			$refpopcol = $colnum;
		}
	}
	
	my %inlinkage;
	open (my $map_handle, "$outdir/$pheno.chr$chr.map") or die "Cannot read $outdir/$pheno.chr$chr.map: $?.\n";
	while ( <$map_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $rsid, $haldane) = split("\t", $_);
		$inlinkage{$rsid} = 1;
	}
	close $map_handle;
	
	open (my $merlinfrq_handle, ">", "$outdir/$pheno.chr$chr.freq") or die "Cannot write to $outdir/$pheno.chr$chr.freq: $?.\n";
	open (my $input_handle, "$refdatadir/freqs/grid$chr.freqs") or die "Cannot read $refdatadir/freqs/grid$chr.freqs: $?.\n";
	while ( <$input_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($SNPid, $rsID, $b37, $ref, $alt, @popfreqs) = split("\t", $_);  	# order of population allele freqs: $AFR, $AMR, $ASN, $EUR, $OVERALL 
		if (exists $inlinkage{$rsID}) {
			my $reffreq = 1 - $popfreqs[$refpopcol-5];
			print $merlinfrq_handle "M $rsID\nA $ref $reffreq\nA $alt $popfreqs[$refpopcol-5]\n";
		}
		
	}
	close $input_handle;
	close $merlinfrq_handle;
}

# perl ~/bin/linkage_scripts/plink2merlin.pl --rawgenodir /net/grc/vol1/mendelian_projects/pastor_uwcmg_et_1/sample_qc/PLINK_100413_0958 --chip ExomeChip --refpop EUR --pheno et --familyedits et.pedchanges.txt --model dominant --outdir /net/grc/vol1/mendelian_projects/pastor_uwcmg_et_1/ngs_analysis/jessica/linkage/merlin_input --interimdir interim_files




################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


plink2merlin.pl - Given the PLINK-formatted output by GenomeStudio, create Merlin-format linkage files.


=head1 SYNOPSIS


perl B<plink2merlin.pl> I<[options]>


=head1 OPTIONS


=over 4

=item B<--rawgenodir> F<directory>	

	path to genotypes in PLINK format as outputted by UWCMG pipeline (likely under project/sample_qc/PLINK*)

=item B<--chip> I<Cytochip|ExomeChip>

	name of genotyping chip
	
=item B<--refpop> I<EUR|AFR|AMR|ASN>

	abbreviation for 1000 Genomes population from which to derive reference allele frequencies

=item B<--pheno> I<name of phenotype>

	name of phenotype for this project (files will use this as a prefix)

=item B<--model> I<dominant|recessive>

	generate a template Merlin .model file for linkage analysis (edit manually to customize parameters)

=item B<--familyedits> F<filename>

	formatted file with edits to be made manually to pedigree (adding missing ancestors, etc)

=item B<--interimdir> F<directory>	

	path to directory to store intermediate files

=item B<--outdir> F<directory>

	path to output directory for finished Merlin format files

=item B<--help> I<help>

	print documentation

=back


=head1 CAVEATS


Remember to edit pheno.model file used by Merlin to customize model of inheritance, causal allele frequency, and penetrance.


=head1 FILES


This script assumes the following files are present in the current directory.


=over 4

=item F<phenotype.updateFID.txt>

	A tab-delimited, 4 column file used by PLINK's --update-ids option to update the family IDs of the genotyped subjects.
	This is necessary because the UWCMG pipeline doesn't output pedigree information with the raw genotypes.  
	The first two columns are the old familyID and subjectID (as listed in the sample_qc/PLINK_*/.ped file).
	The third and fourth columns are the new familyID and subjectID.
	
	Example:  Changes subjectA from family1 to family5 and subjectB from family2 to family6
	family1	subjectA	family5	subjectA
	family2	subjectB	family6	subjectB

=item F<pheno.updateparents.txt>

	A tab-delimited, 4 column file used by PLINK's --update-parents option to update the parent IDs for the genotyped subjects.
	This is necessary because the UWCMG pipeline doesn't output pedigree information with the raw genotypes.  
	The first two columns are the old familyID and subjectID (as listed in the sample_qc/PLINK_*/.ped file).
	The third and fourth columns are the new familyID and subjectID.
	
	Example: Modify the genotype file to include two grandparents (1 and 2) with missing genotypes and to exclude subject 5:
	family1	#1	0	0	1	0
	family1	#2	0	0	2	0
	family1	3	1	2	1	2
	family1	4	1	2	1	2
	family1	!5	1	2	2	0

=item F<familyedits (provided as argument with --familyedits)>

	An optional tab-delimited, 6 column pedigree-format file that describes changes to be made to the genotype data during conversion. 
	This file contains the desired final pedigree information and uses ! to mark subjectIDs that should be excluded and # to mark
	subjectIDs of phantom individuals to be added in (with genotypes listed as missing).  Assumes unique subjectIDs.
	IMPORTANT: If an individual is not listed in this file, they will be excluded from the linkage analysis files.
	NOTE: Merlin and PLINK assume unaffected=1, affected=2, missing=0 but PLINK files from UWCMG pipeline use missing=-9.
	
	Example: Modifies the genotype file to include two grandparents (1 and 2) with missing genotypes and exclude subject 5 completely:
	family1	#1	0	0	1	0
	family1	#2	0	0	2	0
	family1	3	1	2	1	2
	family1	4	1	2	1	2
	family1	!5	1	2	2	0

=back


=head1 EXAMPLES


	perl plink2merlin.pl 
		--rawgenodir /net/grc/vol1/mendelian_projects/pheno/sample_qc/PLINK_100413_0958 
		--chip ExomeChip
		--refpop EUR
		--pheno pheno
		--model dominant
		--familyedits pheno.pedchanges.txt
		--interimdir /net/grc/vol1/mendelian_projects/myphenotype/ngs_analysis/linkage/temp
		--outdir /net/grc/vol1/mendelian_projects/myphenotype/ngs_analysis/linkage/merlin


=head1 AUTHOR


Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)


=cut
