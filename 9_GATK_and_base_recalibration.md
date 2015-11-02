# GATK and Base Recalibration

Or you can go back to realigning indels [here](https://github.com/evansbenj/BIO720/blob/master/8_GATK_realigning_indels.md).

Now that we have the indels realigned, we can proceed with base recalibration.

There are several `GATK` functions that we need to execute to accomplish this. We can do these in one perl script that does the following:
1. Use UnifiedGenotyper to output a **vcf file** with variable sites
2. Use BaseRecalibrator to output a **base recalibration table** based on these sites
3. Use PrintReads to output a **new concatenated bam file** with Recalibrated quality scores
4. Use UnifiedGenotyper to recall the genotypes using the new quality scores, output a new **vcf file**.

``` perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted_realigned.bam files, and 
# make and execute several GATK commandlines with these files that will
# (1) Use UnifiedGenotyper to output a vcf file with variable sites
# (2) Use BaseRecalibrator to output a base recalibration table based on these sites
# (3) Use PrintReads to output a new concatenated bam file with Recalibrated quality scores
# (4) Use UnifiedGenotyper to recall the genotypes using the new quality scores.

my $path_to_reference_genome="~/my_monkey_chromosome/";
my $reference_genome="chrXXX.fa";
my $status;
my @files;

# read in the files within this directory a specific ending
@files = glob("*_sorted_realigned.bam");

# construct a commandline for the UnifiedGenotyper; output a vcf file with only variable sites
my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome." ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline." -out_mode EMIT_VARIANTS_ONLY -o ./nonrecal_varonly_unifiedgenotyper.vcf";

# Execute this command line to call genotypes
$status = system($commandline);

# Now make a new commandline for the BaseRecalibrator, output a table for base recalibration 
$commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R ".$path_to_reference_genome.$reference_genome." ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline." -knownSites ./nonrecal_varonly_unifiedgenotyper.vcf -o recal_data.table";

# Execute this command line
$status = system($commandline);


# Now use PrintReads to output a new concatenated bam file with recalibrated quality scores
my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T PrintReads -R ".$path_to_reference_genome.$reference_genome." ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-BQSR recal_data.table -o concatentated_and_recalibrated_round1.bam";


$status = system($commandline);


# Now use UnifiedGenotyper to recall bases with recalibrated quality scores; output a vcf file 
my $commandline = "java -Xmx3G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ".$path_to_reference_genome.$reference_genome." ";
$commandline = $commandline." -I concatentated_and_recalibrated_round1.bam";
$commandline = $commandline." -o concatentated_and_recalibrated_round1.vcf";


$status = system($commandline);

```
