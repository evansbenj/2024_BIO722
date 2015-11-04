# GATK and Base Recalibration

## Or you can go back to realigning indels [here](https://github.com/evansbenj/BIO720/blob/master/8_GATK_realigning_indels.md).

## Base recalibration and recalling genotypes

Check out a description of the problem [here](http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr).

Now that we have the indels realigned, we can proceed with base recalibration.

There are several `GATK` functions that we need to execute to accomplish this. We can do these in one perl script that does the following:
* Use **UnifiedGenotyper** to output a **vcf file** with variable sites.
* Use **BaseRecalibrator** to output a **base recalibration table** based on these sites.
* Use **PrintReads** to output a **new concatenated bam file** with Recalibrated quality scores.
* Use **UnifiedGenotyper** to call the genotypes using the new quality scores, and output a new **vcf file**.

I wrote another perl script that can do this for us below. In class, please run this on the subset realigned bam files. You could name it `Step_3_base_recalibration.pl`

``` perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted_realigned.bam files, and 
# make and execute several GATK commandlines with these files that will
# (1) Use UnifiedGenotyper to output a vcf file with variable sites
# (2) Use BaseRecalibrator to output a base recalibration table based on these sites
# (3) Use PrintReads to output a new concatenated bam file with Recalibrated quality scores
# (4) Use UnifiedGenotyper to recall the genotypes using the new quality scores; output a new vcf file.

my $path_to_reference_genome="~/my_monkey_chromosome/";
my $reference_genome="chrXXX.fa";
my $status;
my @files;

# read in the files within this directory a specific ending
@files = glob("*_sorted_realigned.bam");

# construct a commandline for the UnifiedGenotyper; output a vcf file with only variable sites
my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}
$commandline = $commandline."-out_mode EMIT_VARIANTS_ONLY -o ./nonrecal_varonly.vcf";

# Execute this command line to call genotypes
$status = system($commandline);

# Now make a new commandline for the BaseRecalibrator, output a table for base recalibration 
$commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R ".$path_to_reference_genome.$reference_genome;
foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-knownSites ./nonrecal_varonly.vcf -o recal_data.table";

# Execute this command line
$status = system($commandline);


# Now use PrintReads to output a new concatenated bam file with recalibrated quality scores
$commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T PrintReads -R ".$path_to_reference_genome.$reference_genome;

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-BQSR recal_data.table -o concatentated_and_recalibrated_round1.bam";


$status = system($commandline);


# Now use UnifiedGenotyper to recall bases with recalibrated quality scores; output a vcf file 
$commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome;
$commandline = $commandline." -I concatentated_and_recalibrated_round1.bam";
$commandline = $commandline." -o recalibrated_round1.vcf";


$status = system($commandline);

```

As previously, you will need to modify this script to point to the path and filename of your chromosome.

## Comparing genotype calls before and after base recalibration

How can we determine whether base recalibration made a difference?  To accomplish this, we can use software called `vcftools` to compare the base calls made before (`nonrecal_varonly.vcf`) and after (`concatentated_and_recalibrated_round1.vcf`) base recalibration.

To use `vcftools`, we need to gzip and index the vcf files we will compare. Please type this:

`/usr/local/tabix/bgzip nonrecal_varonly.vcf`

`/usr/local/tabix/bgzip recalibrated_round1.vcf`

This makes two files with a suffix `gz`.  Now please type this:

`/usr/local/tabix/tabix -p vcf nonrecal_varonly.vcf.gz`

`/usr/local/tabix/tabix -p vcf recalibrated_round1.vcf.gz`

In order for `vcftools` to work we need to make sure it knows where to find the `tabix` program. Please type this:

`export PATH=${PATH}:/usr/local/tabix/`

Now we can use the `vcf-compare` module of `vcftools` to compare these vcf files like this:

`/usr/local/vcftools/src/perl/vcf-compare xxx.vcf.gz yyy.vcf.gz > compare.out`

so you can type this:

`/usr/local/vcftools/src/perl/vcf-compare nonrecal_varonly.vcf.gz recalibrated_round1.vcf.gz > compare.out`

Now check out the output (`more compare.out`), which is basically a Venn diagram of SNPs in each vcf file. Is there a big difference between these files and if so what is the nature of this difference?


# Problem 8. At home, please run this script on your bam files with the complete data.

# OK now we can move on to using `GATK` to filter our data [here](https://github.com/evansbenj/BIO720/blob/master/10_Using_GATK_to_filter_data.md).
