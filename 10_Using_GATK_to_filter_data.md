# Using GATK to filter data

Or go back to GATK and base recalibration [here](https://github.com/evansbenj/BIO720/blob/master/9_GATK_and_base_recalibration.md).

It is often the case, despite our efforts to generate high quality genotype calls, that some genotypes just don't make sense.  For example, we might observe a heterozygous genotype on the male specific portion of the Y chromosome or we might see some genotypes from a female on the Y chromosome. We can easily identify and screen out these sites (i.e. filter them) using `GATK`.

For the purposes of this class, we will first generate a vcf file that has all of the called sites in it (previous vcf files had only variable sites). We will then filter these data based on whether they are near an insertion/deletion (indel), and output a filtered vcf file that has these sites removed. 

Ben will discuss other filtering steps we can do that are relevant to datasets mapped to a complete genome (as opposed to our data which are mapped to individual chromosomes).

Here's an example of a script that can do these steps.  Please copy this and make an executeable file and execute it on the subset data we have been working with. You could name it `Step_4_flag_and_filter.pl`.


``` perl


#!/usr/bin/perl
use warnings;
use strict;

# This script will do the following:
# (1) it will use UnifiedGenotyper to recall bases with recalibrated quality scores
# and output a vcf file with all sites, including variant and homozygous calls.
# (2) It will then use SelectVariants to make another vcf file with only indels in it. 
# (3) then it will use VariantFiltration to mark indels and other low quality sites
# (4) it will use SelectVariants to output a filtered vcf file.



my $path_to_reference_genome="~/my_monkey_chromosome/";
my $reference_genome="chrXXX.fa";
my $status;

# output all sites with UnifiedGenotyper
my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." -I concatentated_and_recalibrated_round1.bam";
$commandline = $commandline." -out_mode EMIT_ALL_CONFIDENT_SITES -o recalibrated_round1_allsites.vcf";
$status = system($commandline);

# make a file with only indels using SelectVariants
$commandline = "java -Xmx2G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T SelectVariants -R .$path_to_reference_genome.$reference_genome; 
$commandline = $commandline." --variant recalibrated_round1_allsites.vcf -selectType INDEL -o indels_only.vcf";
$status = system($commandline);

# flag the vcf file using the indel file
$commandline = "java -Xmx3G -jar GenomeAnalysisTK.jar -T VariantFiltration -R ".$path_to_reference_genome.$reference_genome; 
$commandline = $commandline."-o flagged.vcf --variant recalibrated_round1_allsites.vcf "
$commandline = $commandline." --mask indels_only.vcf --maskName INDEL --maskExtension 10";
$status = system($commandline);

# output a new filtered genotype file using SelectVariants
java -Xmx2g -jar GenomeAnalysisTK.jar -T SelectVariants -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." --variant flagged.vcf -o filtered.vcf -select \'vc.isNotFiltered()\'";
$status = system($commandline);

```

Other examples one could add to the VariantFilteration command line for Ben to discuss:
```
--filterExpression "DP < 5" --filterName "LowCoverage" 
--filterExpression "CHROM == 'chrY' && vc.getGenotype('PF515').isHom()" --filterName "Y_chrom_homoz_filter_for_PF515"
```

# Problem 9.

At home please execute thus script on the complete data.
