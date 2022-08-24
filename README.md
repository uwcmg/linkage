# PLINK2MERLIN(1)       User Contributed Perl Documentation      PLINK2MERLIN(1)



## NAME
       plink2merlin.pl - Given the PLINK-formatted output by GenomeStudio,
       create Merlin-format linkage files.

## SYNOPSIS
       perl plink2merlin.pl [options]

## ARGUMENTS
       --rawgenodir directory
                   path to genotypes in PLINK format as outputted by UWCMG pipeline (likely under project/sample_qc/PLINK*)

       --chip CytoChip|ExomeChip|Omni2_5|OmniExpress
                   name of genotyping chip

       --refpop EUR|AFR|AMR|ASN
                   abbreviation for 1000 Genomes population from which to derive reference allele frequencies

       --pheno name of phenotype
                   name of phenotype for this project (files will use this as a prefix)

       --model dominant|recessive
                   generate a template Merlin .model file for linkage analysis (edit manually to customize parameters)

       --familyedits filename
                   formatted file with final edits to be made manually to pedigree (adding missing ancestors, change phenotype)

       --interimdir directory
                   path to directory to store intermediate files

       --outdir directory
                   path to output directory for finished Merlin format files

       --help help
                   print documentation

## CAVEATS
       Remember to edit pheno.model file used by Merlin to customize model of
       inheritance, causal allele frequency, and penetrance.

## DESCRIPTION
       This script uses PLINK extensively to update family IDs, followed by
       updating parent IDs, followed by excluding SNPs with low call rate,
       zeroing out Mendelian errors, and extracting the grid markers.

## FILES
       This script assumes the following files are present in the current
       directory.

       phenotype.updateFID.txt
                   A tab-delimited, 4 column file used by PLINK's --update-ids option to update the family IDs of the genotyped subjects.
                   This is necessary because the UWCMG pipeline doesn't output pedigree information with the raw genotypes.
                   The first two columns are the old familyID and subjectID (as listed in the sample_qc/PLINK_*/.ped file).
                   The third and fourth columns are the new familyID and subjectID.

                   Example:  Changes subjectA and subjectB from family1 to family5 and subjectC from family2 to family6
                   family1 subjectA        family5 subjectA
                   family1 subjectB        family5 subjectB
                   family2 subjectC        family6 subjectC

       pheno.updateparents.txt
                   A tab-delimited, 4 column file used by PLINK's --update-parents option to update the parent IDs for the genotyped subjects.
                   This is necessary because the UWCMG pipeline doesn't output pedigree information with the raw genotypes.
                   The first two columns are the updated familyID and subjectID as listed in columns 3 and 4 of the phenotype.updateFID.txt file
                   The third and fourth columns are the new paternal ID and maternalID.

                   family5 subjectA        dad     mom
                   family5 subjectB        dad mom
                   family5 dad     0       0
                   family5 mom     0       0
                   family6 subjectC        0       0

       familyedits (provided as argument with --familyedits)
                   An tab-delimited, 6 column PLINK PED format file that describes changes to be made to the pedigree data during conversion.
                   This file contains the desired final pedigree information and uses ! to mark subjectIDs that should be excluded and # to mark
                   subjectIDs of phantom individuals to be added in (with genotypes listed as missing).  Assumes unique subjectIDs.
                   IMPORTANT: If an individual is not listed in this file, they will be excluded from the linkage analysis files.
                   NOTE: Merlin and PLINK assume unaffected=1, affected=2, missing=0 but PLINK files from UWCMG pipeline use missing=-9.
                           Remember to change the affected status in the last column to update phenotype information.

                   Example: Modifies the genotype file to include two grandparents (1 and 2) with missing genotypes and exclude subject 5 completely:
                   family1 #1      0       0       1       0
                   family1 #2      0       0       2       0
                   family1 3       1       2       1       2
                   family1 4       1       2       1       2
                   family1 !5      1       2       2       0

## EXAMPLES
               perl plink2merlin.pl
                       --rawgenodir /net/grc/vol1/mendelian_projects/pheno/sample_qc/PLINK_100413_0958
                       --chip ExomeChip
                       --refpop EUR
                       --pheno pheno
                       --model dominant
                       --familyedits pheno.pedchanges.txt
                       --interimdir /net/grc/vol1/mendelian_projects/myphenotype/ngs_analysis/linkage/temp
                       --outdir /net/grc/vol1/mendelian_projects/myphenotype/ngs_analysis/linkage/merlin

## AUTHOR
       Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)



# perl v5.14.2                      2013-09-20                   PLINK2MERLIN(1)









# MERLIN2BED(1)         User Contributed Perl Documentation        MERLIN2BED(1)



## NAME
       merlin2bed.pl - Given results from Merlin (using --markerNames --tabulate), output a BED file with coordinates of regions meeting a LOD score cutoff and a BED-like file (with extra columns) with more information on the regions.

## SYNOPSIS
       perl merlin2bed.pl [options]

## ARGUMENTS
       --cm2bp *merlin.cm2bp.map
                   PLINK map format to provide the physical positions of variants used in linkage analysis

       --merlintbl .tbl merlin file
                   path to merlin results file (produced when you use --tabulate --markerNames)

       --allchr TRUE|FALSE
                   should this script look across all chromosomes (1:22) with same filename pattern as the --merlintbl argument or just use the one input file

       --cutofflod number
                   minimum LOD score for a marker to be considered part of a linkage peak (peaks will include all markers with LOD>=cutoff plus flanking markers)

       --usechr TRUE|FALSE
                   use the "chr22" notation instead of "22" in the output BED file

       --outprefix filename
                   prefix for files produced by this script

       --help help
                   print documentation

## NOTES
       If plink2merlin.pl was used to generate the input files for Merlin,
       then there should be a file named <phenotype>.merlin.cm2bp.map in the
       same directory as the Merlin files and this file can be used for the
       --cm2bp argument.

       If you have the nextgen data in VCF or BED format, you can Bedtools to
       get the variants within the regions in the output from this script
       using a command like:

       intersectBed -a myphenotype.nextgenvariants.tsv -b linkageresults.bed
       -header > myphenotype.overlaplinkage.tsv

       intersectBed -a myphenotype.exome.vcf -b linkageresults.bed -header >
       myphenotype.overlaplinkage.tsv

## EXAMPLES
       perl merlin2bed.pl
               --cm2bp myphenotype.merlin.map
               --merlintbl myphenotype.dominant.chr22-parametric.tbl
               --cutofflod 1
               --usechr T
               --outprefix temp
               --allchr T

## AUTHOR
       Jessica Chong 


# perl v5.14.2                      2013-09-11                     MERLIN2BED(1)




# GET_HIGHHET_VARIANTS(1)         User Contributed Perl Documentation         GET_HIGHHET_VARIANTS(1)



## NAME
       get_highhet_variants.pl - Given freq and map files (Liz’s format) for a
       given chip, generate a list of markers to be used for a linkage grid
       panel.

## SYNOPSIS
       perl get_highhet_variants.pl [options]

## ARGUMENTS
       --freqdir input file
                   path to directory containing chr$chr.$chip.freq files with population allele frequencies

       --mapdir input file
                   path to directory containing chr$chr.$chip.map files with haldane distsances

       --outdir input file
                   directory to save files produced by this script

       --minfreq input file
                   required minimum allele frequency in any population (if lower, do not use marker in grid)

       --maxfreq input file
                   required maximum allele frequency in any population (if higher, do not use marker in grid)

       --cMdist input file
                   minimum haldane distance (cM) between markers in grid

       --chip ExomeChip|CytoChip
                   name of genotyping chip

       --help help
                   print documentation

## DESCRIPTION
       For each variant in the freq files, check if the alt allele frequency
       is between (minfreq) and (maxfreq) for all populations in the file, if
       the variant passes this check, make sure it has a sex-averaged Haldane
       map position available, then check if marker is at least (cMdist)
       centimorgans away from the previous marker.  If so, add marker to the
       linkage grid list.  If not, continue onto next marker.

## EXAMPLES
               perl get_highhet_variants.pl
                       --freqdir ../../CytoChipComplete/freqs/
                       --mapdir ../../CytoChipComplete/maps/
                       --outdir .
                       --minfreq 0.3
                       --maxfreq 0.7
                       --cMdist 0.5
                       --chip CytoChip

## AUTHOR
       Jessica Chong (jxchong@uw.edu, jxchong@gmail.com)



# perl v5.14.2                      2013-09-20           GET_HIGHHET_VARIANTS(1)
