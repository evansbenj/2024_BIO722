# The Genome Analysis Toolkit (GATK): base recalibration and SNP calling.

## Or you can go back to using `Stacks` to generate files for Phylip and Structure [here](https://github.com/evansbenj/BIO720/blob/master/7_Stacks_and_Structure.md).

So far we have mapped Illumina reads to a reference genome using the programs `samtools` and `bwa` and we have expedited this process with a bash script. We have also used the program `Stacks` to perform various tasks including demultiplexing our data (using the `process_RADtags` program), calculate various polymorphism statistics (using the `populations` program), and output input files for other software such as `Phylip` and `Structure`. 

There are still a few concerns with these approaches and one might view them as `exploratory` analyses for that reason: 
* Once concern is that `Stacks` uses a preset error rate to make SNP calls as opposed to using the quality scores from the fastq files. This is potentially less a less accurate approach than using the fastq quality scores. 
* And even if we did use the fastq quality scores, another concern is that there can be variation among reads in the accuracy of the quality scores. You may recall that fastq files report quality scores on a Phred-scale as described [here](https://en.wikipedia.org/wiki/Phred_quality_score). It turns out that a Phred scaled quality score of "30" actually is not associated with the same confidence level of base call confidence across all sites, and the variation among sites in quality scores is influenced by the genomic context in which the base occurs (for example what bases are next to it) and also the stage of the run that the data is collected (whether early in the datacollection or late). 
* Yet another concern is that each individual in our study were aligned to the reference genome independently. Consequently, if there are insertion/deletion differences between the reference sequence and our study individuals, there may be differences among individuals in the way `bwa` reconciled the mapping of the data to the reference genome, which could generate false SNPs.

The Genome Analysis Toolkit [GATK](https://www.broadinstitute.org/gatk/) provides solutions to these concerns and in general may provide better genotype calls than `Stacks`. Similar to `Stacks`, `GATK` is a suite of programs or functions.  `GATK` is written in java. In class we will use `GATK` to perform base recalibration on our data and to perform genotype calls using the adjusted quality scores concurrently across all of the individuals in our study.


