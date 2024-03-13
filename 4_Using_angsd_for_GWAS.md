# Using angsd for a GWAS

(Or you can go back to using a bash script [here](https://github.com/evansbenj/BIO720/blob/master/3_Lecture_3_Automating_alignment_with_bash.md)).

## Angsd

[Angsd](https://www.popgen.dk/angsd/index.php/ANGSD) is a software suite that can do several analyses, including Genome-wide association studies (GWAS). The goal of a GWAS study is to identify genomic regions where molecular variation (usually single nucleotide polymorphisms, SNPs) are associated with some phenotype of interest. You can learn more about this approach [here](https://en.wikipedia.org/wiki/Genome-wide_association_study).

This software can use bam files as input. So let's make symbolic links to lots of bam files that I made earlier. Please enter your `bam_files` directory and type this:
```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/bams_mapped_to_XLv10_concatscaf_readgroups/* .
```


## Problem 5

With this information in hand, ...


Now lets quickly look into how to generate genotypes from our bam files by clicking [here](https://github.com/evansbenj/2024_BIO722/blob/master/5_genotyping_with_samtools_and_bcftools.md)



