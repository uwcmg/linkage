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
	'familyedits=s' => \$familyedits, 
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
} elsif (!defined $familyedits) {
	pod2usage(-exitval=>2, -verbose=>1, -message => "$0: --familyedits not defined\n");
}


$rawgenodir = remove_trailing_slash_dir($rawgenodir);
$outdir = remove_trailing_slash_dir($outdir);
$interimdir = remove_trailing_slash_dir($interimdir);

my ($updatecMfile, $updatecMsnplist, $refdatadir);


print "\n";
print "Using $genotypechip markers and frequencies because you specified genotypechip=$genotypechip\n";	
if ($genotypechip =~ /ExomeChip/i) {
	$refdatadir = "$mendeliandir/ExomeChipLinkage";
	$updatecMfile = 'allchr.ExomeChip.updatecM.nodups.txt';
	$updatecMsnplist = 'allchr.ExomeChip.updatecM.snplist';
	if (system("cut -f1,7 $mendeliandir/ExomeChipCompleteMaps/chr*.ExomeChip.map > $interimdir/allchr.ExomeChip.updatecM.txt") != 0) {
		die "Can't extract haldane values from $mendeliandir/ExomeChipCompleteMaps/chr*.ExomeChip.map: $?";
	}
	if (system("cut -f1,7 $mendeliandir/ExomeChipLinkage/maps/grid*.map >> $interimdir/allchr.ExomeChip.updatecM.txt") != 0) {
		die "Can't extract haldane values from $mendeliandir/ExomeChipLinkage/maps/grid*.map: $?";
	}
	if (system("sort -k1 $interimdir/allchr.ExomeChip.updatecM.txt | uniq > $interimdir/allchr.ExomeChip.updatecM.nodups.txt") != 0) {
		die "Can't get only unique SNPs from $interimdir/allchr.ExomeChip.updatecM.txt: $?";
	}
	if (system("cut -f1 $interimdir/allchr.ExomeChip.updatecM.txt > $interimdir/allchr.ExomeChip.updatecM.snplist") != 0) {
		die "Can't extract list of SNPs with haldane values from $interimdir/allchr.ExomeChip.updatecM.snplist: $?";
	}
	if (system("cut -f2 $refdatadir/freqs/grid*.freqs > $interimdir/gridrsIDvariantsonly.snplist")) {
		die "Can't extract rsIDs for linkage analysis with cM values from $refdatadir/freqs/grid*.freqs: $?";
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


# generate list of variants with rsIDs that are in the raw genotype data,
#   cleanup the name to be only the rsID itself, create update_rsIDS.txt
create_update_rsID_names($interimdir, $rawgenodir, $rawgenomap, $refdatadir);


print "... creating PLINK files: extracting SNPs with rsIDs, updating genetic maps, extracting SNPs with call rate >= 95%\n";
# update raw genotype files with rsIDs and create new PLINK files with only those variants
`plink --file $rawgenodir/$rawgenostem --update-map $interimdir/update_rsIDs.txt --update-name --make-bed --out $interimdir/$pheno.gridvar`;
`plink --bfile $interimdir/$pheno.gridvar --extract $interimdir/rsIDvariantsonly.snplist --make-bed --out $interimdir/$pheno.gridrsIDs`;

# update PLINK files with genetic map information
`plink --bfile $interimdir/$pheno.gridrsIDs --update-map $interimdir/$updatecMfile --update-cm --extract $interimdir/$updatecMsnplist --make-bed --out $interimdir/$pheno.haldane`;
### verify all variants have genetic map updated: "0 in data but not in [ /nfs/home/jxchong/lib/allchr.ExomeChip.updatecM.txt ]"

# do some basic QC: minimum 95% genotyping rate
`plink --bfile $interimdir/$pheno.haldane --geno 0.05 --nonfounders --make-bed --out $interimdir/$pheno.callrate95`;

print "... creating PLINK files: updating family information\n";
# update family ID, phenotype, parental information based on files
`plink --bfile $interimdir/$pheno.callrate95 --update-ids $pheno.updateFID.txt --make-bed --out $interimdir/$pheno.updateFID`;
`plink --bfile $interimdir/$pheno.updateFID --update-parents $pheno.updateparents.txt --make-bed --out $interimdir/$pheno.updateparents`;

print "... creating PLINK files: determining SNPs that need allele flipping\n";
# extract only variants in the grid files
`plink --bfile $interimdir/$pheno.updateparents --extract $interimdir/gridrsIDvariantsonly.snplist --make-bed --out $interimdir/$pheno.gridonly`;

# flip alleles in genotype file to match the 1000 Genomes reference file
makefliplist("$pheno.gridonly.bim", $refdatadir);
`plink --bfile $interimdir/$pheno.gridonly --flip $interimdir/rsIDs.toflip.txt --exclude $interimdir/rsIDs.toexclude.txt --make-bed --out $interimdir/$pheno.flipped`;

# do mendelian error check and zero out all probematic genotypes
`plink --bfile $interimdir/$pheno.flipped --me 1 1 --set-me-missing --missing-phenotype 0 --recode --out $interimdir/$pheno.updateparents.me1-1`;
# need to use --transpose so it is easier to sort by genetic position???


copy("$interimdir/$pheno.updateparents.me1-1.map", "$interimdir/$pheno.merlin.map") or die "Failed to copy $interimdir/$pheno.updateparents.me1-1.map to $interimdir/$pheno.merlin.map\n";
copy("$interimdir/$pheno.merlin.map", "$outdir/$pheno.merlin.cm2bp.map") or die "Failed to copy $interimdir/$pheno.merlin.map to $outdir/$pheno.merlin.cm2bp.map\n";


# add in any extra ancestors and exclude people if desired
# if (defined $familyedits) {
	print "\tPedigree edits provided ($familyedits); making changes.\n";
	if (system("perl /net/grc/vol1/mendelian_projects/mendelian_analysis/module_linkage/edit_linkage_family_info.pl --edits $familyedits --pedfile $interimdir/$pheno.updateparents.me1-1.ped --out $interimdir/$pheno.merlin.ped") != 0) {
		die "Could not run: perl /net/grc/vol1/mendelian_projects/mendelian_analysis/module_linkage/edit_linkage_family_info.pl --edits $familyedits --pedfile $interimdir/$pheno.updateparents.me1-1.`ped --out $interimdir/$pheno.merlin.ped: $?";
	}
# } else {
# 	copy("$interimdir/$pheno.updateparents.me1-1.ped", "$interimdir/$pheno.merlin.ped") or die "Failed to copy $interimdir/$pheno.updateparents.me1-1.ped to $interimdir/$pheno.merlin.ped\n";
# }


# create a default model file if one doesn't already exist
if (! -e "$outdir/$pheno.model") {
	if ($model =~ /dominant/i) {
		`echo '$pheno 0.001 0,0.95,1.0 DOES_THIS_MODEL_MAKE_SENSE_(q=0.001,_f=0,0.95,1.0)' | cat > $outdir/$pheno.model`;
	} elsif ($model =~ /recessive/i) {
		`echo '$pheno 0.005 0,0,1.0 DOES_THIS_MODEL_MAKE_SENSE_(q=0.005,_f=0,0,1.0)' | cat > $outdir/$pheno.model`;
	} else {
		`echo '$pheno 0.005 0,0,1.0 DOES_THIS_MODEL_MAKE_SENSE_(q=0.005,_f=0,0,1.0)' | cat > $outdir/$pheno.model`;
	}
}


# Output linkage file per chromosome
print "... making linkage files\n";
for (my $chr=1; $chr<=22; $chr++) {
	print "\tfor chromosome $chr\n";
	# create ped file
	`plink --file $interimdir/$pheno.merlin --chr $chr --missing-phenotype 0 --recode --out $interimdir/$pheno.merlin.chr$chr`;
	move("$interimdir/$pheno.merlin.chr$chr.ped", "$outdir/$pheno.chr$chr.ped") or die "Failed to move $interimdir/$pheno.merlin.chr$chr.ped to $outdir/$pheno.chr$chr.ped\n";
	
	# create map file 
	`cut -f1-3 $interimdir/$pheno.merlin.chr$chr.map > $outdir/$pheno.chr$chr.map`;
	
	# create frequency file
	create_frqfile($chr, $refdatadir, $refpop, $pheno, $outdir);
	# allele flipping to account for different strands

	# create .dat file
	`echo 'A $pheno' > $outdir/$pheno.chr$chr.dat`;
	`perl -ane \'print \"M \$F[1]\\n\";\' $interimdir/$pheno.merlin.chr$chr.map >> $outdir/$pheno.chr$chr.dat`;
}

print "Merlin format files $pheno.chr*.map $pheno.chr*.ped $pheno.chr*.freq $pheno.chr*.dat $pheno.model are ready in $outdir.\n";
print "$pheno.model should be edited to fit your model of inheritance\n\n";



# create job submission script to run linkage
open (my $submit_handle, ">", "$outdir/$pheno.sublinkage.sh") or die "Cannot write to $outdir/$pheno.sublinkage.sh: $?.\n";
print $submit_handle '#$ -S /bin/bash'."\n";
print $submit_handle '#$ -t 1-22'."\n";
print $submit_handle '#$ -o .'."\n";
print $submit_handle '#$ -e .'."\n";
print $submit_handle '#$ -cwd'."\n\n";
print $submit_handle 'module load modules modules-init modules-gs gsits-util/latest'."\n";
print $submit_handle 'module load plink/latest perl/5.14.2 merlin/latest'."\n";
print $submit_handle "\n";
print $submit_handle "merlin -d $pheno.chr\${SGE_TASK_ID}.dat -p $pheno.chr\${SGE_TASK_ID}.ped -m $pheno.chr\${SGE_TASK_ID}.map -f $pheno.chr\${SGE_TASK_ID}.freq --model $pheno.model --pdf --prefix $pheno.$model.chr\${SGE_TASK_ID} --tabulate --markerNames\n\n\n";

print $submit_handle "### in case you need to manually re-run linkage on one of the chromosomes:\n";
for (my $i=1; $i<=22; $i++) {
	print $submit_handle "## merlin -d $pheno.chr$i.dat -p $pheno.chr$i.ped -m $pheno.chr$i.map -f $pheno.chr$i.freq --model $pheno.model --pdf --prefix $pheno.$model.chr$i --tabulate --markerNames\n";
}
close $submit_handle;

print "To run linkage: qsub $outdir/$pheno.sublinkage.sh\n";
print "To create BED file summarizing linkage results, use /net/grc/vol1/mendelian_projects/mendelian_analysis/module_linkage/merlin2bed.pl\n";
print "Linkage plots are created by Merlin: see $outdir/$pheno.$model.chr*.pdf\n";
print "Linkage result tables are created by Merlin: see $outdir/$pheno.$model.chr*-parametric.tbl\n";






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

sub create_update_rsID_names {
	# use the rsIDs in the original exported PLINK map file and also try to match with grid files based on position and alleles
	my ($interimdir, $rawgenodir, $rawgenomap, $refdatadir) = @_;
	
	# make plink bed file so it's easier to access the allele and position information
	`plink --file $rawgenodir/$rawgenostem --make-bed --out $interimdir/originalgenotypes`;
	
	# read through bim file to get alleles and position; get likely rsID match from the reference grid files
	## note that grid files only contain ~5k markers
	my %gridref;
	for (my $chr=1; $chr<=22; $chr++) {
		open (my $grid_ref_handle, "$refdatadir/freqs/grid$chr.freqs") or die "Cannot read $refdatadir/freqs/grid$chr.freqs: $!.\n";
		while ( <$grid_ref_handle> ) {
			$_ =~ s/\s+$//;					# Remove line endings
			my ($SNPid, $rsID, $b37, $ref, $alt, @popfreqs) = split("\t", $_); 
			$gridref{$SNPid} = $rsID;
			# $gridref{"$chr:$b37"}{'ref'} = $ref;
			# $gridref{"$chr:$b37"}{'alt'} = $alt;
			# $gridref{"$chr:$b37"}{'rsID'} = $rsID;
		}
		close $grid_ref_handle;
	}
	
	open (my $update_rsid_handle, ">", "$interimdir/update_rsIDs.txt") or die "Cannot write to $interimdir/update_rsIDs.txt: $!.\n";
	open (my $raw_map_handle, "$interimdir/originalgenotypes.bim") or die "Cannot read $interimdir/originalgenotypes.bim: $!.\n";
	while ( <$raw_map_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $varname, $cM, $bp, @alleles) = split("\t", $_);
		if ($varname =~ 'rs') {
			my $rsid = $varname;
			$rsid =~ s/exm-//;
			$rsid =~ s/newrs/rs/;
			print $update_rsid_handle "$varname\t$rsid\n";
		} elsif (defined $gridref{$varname}) {
			print $update_rsid_handle "$varname\t$gridref{$varname}\n";
		}
	}
	close $raw_map_handle;
	close $update_rsid_handle;
	
	if (system("cut -f2 $interimdir/update_rsIDs.txt > $interimdir/rsIDvariantsonly.snplist")) {
		die "Can't extract all SNPs with rsIDs: $?";
	}  
}

sub makefliplist {
	my ($bimfile, $refdatadir) = @_;	

	my %ref_alleles;
	my ($grid_nvar, $reject_nvar, $ambig_nvar, $flip_nvar, $bim_nvar) = ((0) x 5);
	for (my $chr=1; $chr<=22; $chr++) {
		open (my $ref_allele_handle, "$refdatadir/freqs/grid$chr.freqs") or die "Cannot read $refdatadir/freqs/grid$chr.freqs: $?.\n";
		while (<$ref_allele_handle>) {
			$_ =~ s/\s+$//;					# Remove line endings
			my ($SNPid, $rsID, $b37, $ref, $alt, @popfreqs) = split("\t", $_); 
			$ref_alleles{$rsID}{'ref'} = $ref;
			$ref_alleles{$rsID}{'alt'} = $alt;
			$grid_nvar++;
		}
		close $ref_allele_handle;
	}
	
	open (my $fliplist_handle, ">", "$interimdir/rsIDs.toflip.txt") or die "Cannot write to $interimdir/rsIDs.toflip.txt: $?.\n";
	open (my $ambiguouslist_handle, ">", "$interimdir/rsIDs.toexclude.txt") or die "Cannot write to $interimdir/rsIDs.toexclude.txt: $?.\n";
	open (my $bim_handle, "$interimdir/$bimfile") or die "Cannot read $interimdir/$bimfile: $?.\n";
	while (<$bim_handle>) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $rsID, $cM, $b37, $a1, $a2) = split("\t", $_);
		if ($chr > 22) {
			last;						# only handling autosomal chromosomes for now
		}
		$bim_nvar++;
		my $action = needsflip($a1, $a2, $ref_alleles{$rsID}{'ref'}, $ref_alleles{$rsID}{'alt'});
		if ($action eq 'flip') {
			print $fliplist_handle "$rsID\n";
			$flip_nvar++;
		} elsif ($action eq 'ambiguous' || $action eq 'monomorphic') {
			print $ambiguouslist_handle "$rsID\n";
			$ambig_nvar++;
		} elsif ($action eq 'weird') {
			print STDERR "$rsID on chr$chr is not ambiguous (AT or GC SNP) but doesn't need to be flipped either.\n";
			print $ambiguouslist_handle "$rsID\n";
			$reject_nvar++;
		}
	}
	close $bim_handle;	
	close $ambiguouslist_handle;	
	close $fliplist_handle;
	
	print "... ... out of $bim_nvar variants with genotypes and in the linkage grid files, rejecting $ambig_nvar ambiguous and $reject_nvar weird variants, and flipping $flip_nvar variants\n";
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
	my ($chr, $refdatadir, $refpop, $pheno, $outdir) = @_;
	
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
	open (my $map_handle, "$outdir/$pheno.chr$chr.map") or die "Cannot read 1 $outdir/$pheno.chr$chr.map: $?.\n";
	while ( <$map_handle> ) {
		$_ =~ s/\s+$//;					# Remove line endings
		my ($chr, $rsid, $haldane) = split("\t", $_);
		$inlinkage{$rsid} = 1;
	}
	close $map_handle;
	
	open (my $merlinfrq_handle, ">", "$outdir/$pheno.chr$chr.freq") or die "Cannot write to $outdir/$pheno.chr$chr.freq: $?.\n";
	open (my $input_handle, "$refdatadir/freqs/grid$chr.freqs") or die "Cannot read 2 $refdatadir/freqs/grid$chr.freqs: $?.\n";
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





################################################################################################################
############################################ Documentation #####################################################
################################################################################################################


=head1 NAME


plink2merlin.pl - Given the PLINK-formatted output by GenomeStudio, create Merlin-format linkage files.


=head1 SYNOPSIS


perl B<plink2merlin.pl> I<[options]>


=head1 ARGUMENTS


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

	formatted file with final edits to be made manually to pedigree (adding missing ancestors, change phenotype)

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
	
	Example:  Changes subjectA and subjectB from family1 to family5 and subjectC from family2 to family6
	family1	subjectA	family5	subjectA
	family1	subjectB	family5	subjectB
	family2	subjectC	family6	subjectC

=item F<pheno.updateparents.txt>

	A tab-delimited, 4 column file used by PLINK's --update-parents option to update the parent IDs for the genotyped subjects.
	This is necessary because the UWCMG pipeline doesn't output pedigree information with the raw genotypes.  
	The first two columns are the familyID and subjectID as listed in the phenotype.updateFID.txt file
	The third and fourth columns are the new paternal ID and maternalID.
	
	family5	subjectA	dad	mom
	family5	subjectB	dad mom
	family5	dad	0	0
	family5	mom	0	0
	family6	subjectC	0	0

=item F<familyedits (provided as argument with --familyedits)>

	An tab-delimited, 6 column PLINK PED format file that describes changes to be made to the pedigree data during conversion.
	This file contains the desired final pedigree information and uses ! to mark subjectIDs that should be excluded and # to mark
	subjectIDs of phantom individuals to be added in (with genotypes listed as missing).  Assumes unique subjectIDs.
	IMPORTANT: If an individual is not listed in this file, they will be excluded from the linkage analysis files.
	NOTE: Merlin and PLINK assume unaffected=1, affected=2, missing=0 but PLINK files from UWCMG pipeline use missing=-9.
		Remember to change the affected status in the last column to update phenotype information.
	
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
