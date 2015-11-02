# The Genome Analysis Toolkit (GATK): indel realignment.

## Or you can go back to using `Stacks` to generate files for Phylip and Structure [here](https://github.com/evansbenj/BIO720/blob/master/7_Stacks_and_Structure.md).

So far we have mapped Illumina reads to a reference genome using the programs `samtools` and `bwa` and we have expedited this process with a bash script. We have also used the program `Stacks` to perform various tasks including demultiplexing our data (using the `process_RADtags` program), calculate various polymorphism statistics (using the `populations` program), and output input files for other software such as `Phylip` and `Structure`. 

There are still a few concerns with these approaches and one might view them as `exploratory` analyses for that reason: 
* Once concern is that `Stacks` uses a preset error rate to make SNP calls as opposed to using the quality scores from the fastq files. This is potentially less a less accurate approach than using the fastq quality scores. 
* And even if we did use the fastq quality scores, another concern is that there can be variation among reads in the accuracy of the quality scores. You may recall that fastq files report quality scores on a Phred-scale as described [here](https://en.wikipedia.org/wiki/Phred_quality_score). It turns out that a Phred scaled quality score of "30" actually is not associated with the same confidence level of base call confidence across all sites, and the variation among sites in quality scores is influenced by the genomic context in which the base occurs (for example what bases are next to it) and also the stage of the run that the data is collected (whether early in the datacollection or late). 
* Another concern is that each individual in our study were aligned to the reference genome independently. Consequently, if there are insertion/deletion differences between the reference sequence and our study individuals, there may be differences among individuals in the way `bwa` reconciled the mapping of the data to the reference genome, which could generate false SNPs.
* And even if we did recalibrate, realign, and recall our genotypes, we still might want to filter regions that might be unreliable (for example regions near indels, or regions with weird genotype calls such as females with reads on the Y chromosome).

The Genome Analysis Toolkit [GATK](https://www.broadinstitute.org/gatk/) provides solutions to these concerns and in general may provide better genotype calls than `Stacks`. Similar to `Stacks`, `GATK` is a suite of programs or functions.  `GATK` is written in java. In class we will use `GATK` to perform base recalibration on our data and to perform genotype calls using the adjusted quality scores concurrently across all of the individuals in our study.

# A brief note about Ben's approach to using `GATK` (and doing analyses in general)

An important part of any study is its ability to be replicated by others. In fact, a journal could reject your paper if there are insufficient details about your methodology to allow for someone to repeat the analysis. For this reason it is minimally crucial to write down every step that you perform in your analysis, including the version of the software you use, and it is adviseable to make available the scripts you use to the community. The latter point can easily be done with databases such as [dryad](http://datadryad.org/). For describing in detail the steps of an analysis, I find it useful to have every step performed by a script. For the purposes of this class, we will perform one step with one script (we will use `perl` scripts), but obviously once could do multiple steps in one script. As you will see (and have seen with the bash script we used earlier), this also makes the building of commandlines quite simple and saves you from typing out sample names and directories by hand.

# Realigning insertion/deletion (indels) using `GATK`

Before we begin, let's first check out again how a vcf file displays information about indels [here](http://samtools.github.io/hts-specs/VCFv4.2.pdf).

We will use `GATK` to identify indels that may be associated with inappropriate mapping differences among the individuals in our study, and then realign them across all individuals. This is done in two steps.  The first uses the `GATK` function called `RealignerTarget` to identify indels in our data. This function produces a text file that has information about the locations of all indels in any individual relative to the reference genome. Because the second step takes a few dozen minutes to run with the full dataset, for class we will work with the subset datasets that you previously mapped to your chromosome.

Here is a perl script that will execute the `RealignerTarget` function in `GATK` on our data:

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute a GATK commandline on these files.  

my $path_to_reference_genome="~/my_monkey_chromosome/";
my $reference_genome="chrXXX.fa";
my @files;
my $status;
   
@files = glob("*_sorted.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator ";

foreach(@files){
    $commandline = $commandline." -I ".$_." ";
}

$commandline = $commandline."-R ".$path_to_reference_genome.$reference_genome." -o forIndelRealigner.intervals";

$status = system($commandline);
```

Please copy and paste this script, change the permissions to allow it to be executable, and execute it on your samples from within the directory that contains your sorted bam files. I suggest naming your scripts using a sensible system, such as with descriptions and numbers in ascending order.  For example, you could name this script `Step_1_execute_GATK_RealignerTargetCreator.pl`. You will need to adjust the `$reference_genome` variable to match your chromosome.

In this script the `@files = glob("*_sorted.bam");` command uses the `glob` function to look for all files in your directory with that end with "_sorted.bam" and then it adds them to an array called `@files`. The `$status = system($commandline);` line executes the test stored in the $commandline variable. The commandline executes the java `.jar` file and allocates a maximum of 1 Gb of memory for the Java virtual machine (`-Xmx1G`).

When it is done, please check out the file it made like this:

`more forIndelRealigner.intervals`

You shoudl see a list of intervals coordinates following the chromosome of your reference chromosome.

With the indel text file, we can then use a function called `IndelRealigner`, which takes as input this `vcf` file to realign bases when possible an minimize mis-called SNPs.

Here is a perl script that executes the `IndelRealigner` function:

```perl
#!/usr/bin/perl                                                                                                                      
use warnings;
use strict;

# This script will read in the *_sorted.bam file names in a directory, and                                                           
# make and execute a GATK commandline on these files.                                                                                

my $path_to_reference_genome="~/my_monkey_chromosome/";
my $reference_genome="chrXXX.fa";
my $status;
my @files;

@files = glob("*_sorted.bam");

my $commandline = "java -Xmx1G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T IndelRealigner ";

foreach(@files){
    $commandline = $commandline." -I ".$_;
}

$commandline = $commandline." -R ".$path_to_reference_genome.$reference_genome." --targetIntervals forIndelRealigner.intervals --nWayOut _realigned.bam";

$status = system($commandline);

```

As above, please copy and paste this script, make it executable, and execute it. You could name this script `Step_2_execute_GATK_IndelRealigner.pl`. Please don't forget to modify the name of your reference chromosome as appropriate.

If your run completed successfully, you should see new bam files that have the ending `*_sorted_realigned.bam`. These files should be about the same size as the `*sorted.bam` files.  Please check this by typing `ls -l`.

# Problem 7. OK please do these steps at home on the full data set.

# OK, now let's move on to base recalibration with GATK [here](https://github.com/evansbenj/BIO720/blob/master/9_GATK_and_base_recalibration.md).
