#!perl
#
# Description:
#
# Usage: perl plink2merlin.pl <PLINK PED/MAP PREFIX> <OUPUT PREFIX>
#
#
# Created by Jessica on 2010-05-13

use warnings;

if (@ARGV < 2) {
	print "\nUsage: perl plink2merlin_bychr.pl <DISEASE NAME NO SPACES> <just \"convert\" to merlin format or run \"merlin\" or do \"both\" or \"print\" commands to run merlin or run only \"parametric / nonparametric\"> TARGETFOLDER\n\n";
	# print qq(\tplink --file hutt500k.masel05 --keep *.findivs --pheno *.findivs --recode12 --out *.premerlin)."\n\tmv *premerlin* ../linkage\n";
	exit;
}

my $pheno = $ARGV[0];
my $action = $ARGV[1];
my $targetfolder = $ARGV[2];
$targetfolder =~ s/\///;
my $modelname;
if (exists $ARGV[3]) {
	$modelname = $ARGV[3];
}

if ($action eq 'convert' or $action eq 'both') {
	print "... making plink file for $pheno\n";
	
	# my $prefix = "$pheno.premerlin";
	# 
	# Make plink file for the phenotype
	if (! -e "/cchome/jxchong/plink_format/$pheno.premerlin.ped") {
		`/apps/rh5_64/plink --file /clusta/jxchong/plink_format/huttqc.masel05 --keep ~/plink_format/$pheno.findivs --pheno ~/plink_format/$pheno.findivs --recode12 --out ~/plink_format/$pheno.premerlin`;
	}
	
	print "... making plink files\n";
	# Output plink file per chromosome
	for (my $i=1; $i<=22; $i++) {
		print "\tfor chromosome $i\n";
		`/apps/rh5_64/plink --file ~/plink_format/$pheno.premerlin --remove ~/plink_format/AffyFindivs_05-11_temporary_exclusions.txt --chr $i --recode --out ~/plink_format/$pheno.premerlin.chr$i`;
		`mv ~/plink_format/$pheno.premerlin.chr$i.* ~/linkage`;
		# Rewrite family information in plink files (assumes this is necessary)
	 		# `/usr/bin/perl edit_linkage_family_info.pl $targetfolder/$pheno.linkagefamily.txt $pheno.premerlin.chr$i.ped > $pheno.merlin.chr$i.ped`;
			`cp $pheno.premerlin.chr$i.ped $pheno.merlin.chr$i.ped`;
		# `cp $pheno.premerlin.chr$i.ped $pheno.original.chr$i.ped`; ## DEBUG
		`rm $pheno.premerlin.chr$i.ped`;
		`sort $pheno.merlin.chr$i.ped > $pheno.premerlin.chr$i.ped`;
		`rm $pheno.merlin.chr$i.ped`;
	}
	
	print "... converting plink files to merlin format\n";
	
	for (my $i=1; $i<=22; $i++) {
		print "\tfor chromosome $i\n";
		`/usr/bin/perl -ane \'print \"\$F[0]\\t\$F[1]\\t\".(\$F[3]\/1000000).\"\\n\";\' $pheno.premerlin.chr$i.map > $pheno.chr$i.map`;
	
		`/usr/bin/perl add_cM_pos.pl $pheno.chr$i.map > $pheno.chr$i.map2`;
	
		# Frequency file generated for each chromosome
		`/usr/bin/perl -ane \'\$major = (1-\$F[4]); \$minor = \$F[4]; print "M \$F[1]\\nF \$minor\\nF \$major\\n\";\' ~/linkage/freqfiles/huttqc.masel05.chr$i.frq > $pheno.chr$i.freq`;
	
		`echo 'A $pheno' | cat > $pheno.chr$i.dat`;
	
		`/usr/bin/perl -ane \'print \"M \$F[1]\\n\";\' $pheno.premerlin.chr$i.map >> $pheno.chr$i.dat`;
	
		`rm $pheno.premerlin.chr$i.map`;
		`rm $pheno.premerlin.chr$i.log`;
		`mv $pheno.chr$i.map2 $pheno.chr$i.map`;
		`mv $pheno.premerlin.chr$i.ped $pheno.chr$i.ped`;
	}
	
	`rm $pheno.premerlin.chr*.nof`;
		
	`echo '$pheno 0.05 0,0,1.0 Recessive_Model' | cat > $pheno.model`;
	
	`mv $pheno.chr* $targetfolder/`;
	`mv $pheno.model $targetfolder/`;
	
	print "Merlin format files $pheno.chr*.map $pheno.chr*.ped $pheno.chr*.freq $pheno.chr*.dat $pheno.model are ready.\n\n";
	
	if ($modelname) {
		open (ALL, ">$targetfolder/allchr.$modelname.sh") or die "Cannot open $targetfolder/allchr.$modelname.sh file.\n";
		print ALL "#!/bin/bash\n\n";
		print "Making qsub .sh file for chromosome: ";
		for (my $i=1; $i<=22; $i++) {
			print "$i..";
			open (FILE, ">$targetfolder/chr$i.$modelname.sh") or die "Cannot open $targetfolder/chr$i.$modelname.sh file.\n";
				print FILE "#!/bin/sh\n#\$ -cwd\n#\$ -S /bin/bash\n\n";
				print FILE qq(/apps/rh5_64/merlin -d $pheno.chr$i.dat -p $pheno.chr$i.ped -m $pheno.chr$i.map -f $pheno.chr$i.freq --tabulate --model $pheno.$modelname.model --markerNames --pdf --prefix $pheno.$modelname.chr$i.pl)."\n";
			close FILE;
			print ALL "/opt/gridengine/bin/lx26-amd64/qsub -V -cwd chr$i.$modelname.sh\n";
		}
		close ALL;
		`chmod u+x $targetfolder/allchr.$modelname.sh`;
		print "\nRun: ./$targetfolder/allchr.$modelname.sh\n";
		print "If not using a recessive model, don't forget to create file $pheno.$modelname.model\n\n";
	} else {
		open (ALL, ">$targetfolder/allchr.sh") or die "Cannot open $targetfolder/allchr.sh file.\n";
		print ALL "#!/bin/bash\n\n";
		print "Making qsub .sh file for chromosome: ";
		for (my $i=1; $i<=22; $i++) {
			print "$i..";
			open (FILE, ">$targetfolder/chr$i.sh") or die "Cannot open $targetfolder/chr$i.sh file.\n";
				print FILE "#!/bin/sh\n#\$ -cwd\n#\$ -S /bin/bash\n\n";
				print FILE qq(/apps/rh5_64/merlin -d $pheno.chr$i.dat -p $pheno.chr$i.ped -m $pheno.chr$i.map -f $pheno.chr$i.freq --tabulate --model $pheno.model --markerNames --pdf --prefix $pheno.chr$i.pl)."\n";
			close FILE;
			print ALL "/opt/gridengine/bin/lx26-amd64/qsub -V -cwd chr$i.sh\n";
		}
		`chmod u+x $targetfolder/allchr.sh`;
		print "\nRun: ./$targetfolder/allchr.sh\n";
	}
}


# if ($action eq 'merlin' or $action eq 'both' or $action eq 'parametric') {
# 	print "... running merlin\n";
# 	# system(qq(merlin -d $pheno.dat -p $pheno.ped -m $pheno.map -f $pheno.freq --npl --pairs --exp --tabulate --pdf));
# 	# system(qq(mv merlin.pdf $pheno.npl.pdf));
# 	# system(qq(mv merlin-nonparametric.tbl $pheno.npl-table.txt));
# 
# 	system(qq(merlin -d $pheno.dat -p $pheno.ped -m $pheno.map -f $pheno.freq --tabulate --model $pheno.model --pdf --prefix $pheno.pl));
# 	system(qq(mv merlin.pdf $pheno.pl.pdf));
# 	system(qq(mv merlin-parametric.tbl $pheno.pl-table.txt));
# }
# 
# if ($action eq 'merlin' or $action eq 'both' or $action eq 'nonparametric') {
# 	print "... running merlin\n";
# 	system(qq(merlin -d $pheno.dat -p $pheno.ped -m $pheno.map -f $pheno.freq --npl --exp --tabulate --pdf --prefix $pheno.npl));
# 	system(qq(mv merlin.pdf $pheno.npl.pdf));
# 	system(qq(mv merlin-nonparametric.tbl $pheno.npl-table.txt));
# }

if ($action eq 'print') {
	for (my $i=1; $i<=22; $i++) {
		# print qq(merlin -d $pheno.chr$i.dat -p $pheno.chr$i.ped -m $pheno.chr$i.map -f $pheno.chr$i.freq --npl --exp --tabulate --markerNames --pdf --prefix $pheno.chr$i.npl)."\n";

		# print qq(mv merlin.pdf $pheno.npl-markerNames.pdf)."\n";
	
		# print qq(mv merlin-nonparametric.tbl $pheno.npl-markerNames-table.txt)."\n";

		print qq(/apps/rh5_64/merlin -d $pheno.chr$i.dat -p $pheno.chr$i.ped -m $pheno.chr$i.map -f $pheno.chr$i.freq --tabulate --model $pheno.model --markerNames --pdf --prefix $pheno.chr$i.pl)."\n";

		# print qq(mv merlin.pdf $pheno.pl-markerNames.pdf)."\n";

		# print qq(mv merlin-parametric.tbl $pheno.pl-markerNames-table.txt)."\n";
	}
}
