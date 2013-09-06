#!perl
#
# Description:
#
# Usage: perl untitled.pl
#
#
# Created by Jessica on 2011-11-21

use strict;
use warnings;


print "Making qsub .sh file for chromosome: ";
for (my $i=1; $i<=22; $i++) {
	# print "$i..";
	open (FILE, ">/clusta/jxchong/ghayda/linkage/merlinfiles/chr$i.sh") or die "Cannot open /clusta/jxchong/ghayda/linkage/merlinfiles/chr$i.CLBA.sh file.\n";
		print FILE "#!/bin/bash\n#\$ -cwd\n#\$ \n\n";
		print FILE qq(/apps/rh5_64/merlin -d CLBA.final.30k.merlin.chr$i.dat -p CLBA.final.30k.merlin.chr$i.ped -m CLBA.final.30k.merlin.chr$i.map --model CLBA.model --tabulate --markerNames --pdf --prefix CLBA.veryrare.chr$i.pl)."\n";
	close FILE;
	print "qsub chr$i.sh\n";
}

