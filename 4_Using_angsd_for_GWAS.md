# Using angsd for a GWAS

(Or you can go back to using a bash script [here](https://github.com/evansbenj/BIO720/blob/master/3_Lecture_3_Automating_alignment_with_bash.md)).

## Angsd

[Angsd](https://www.popgen.dk/angsd/index.php/ANGSD) is a software suite that can do several analyses, including Genome-wide association studies (GWAS). The goal of a GWAS study is to identify genomic regions where molecular variation (usually single nucleotide polymorphisms, SNPs) are associated with some phenotype of interest. You can learn more about this approach [here](https://en.wikipedia.org/wiki/Genome-wide_association_study).

This software can use bam files as input. So let's make symbolic links to lots of bam files that I made earlier. Please make sure you are in your `bam_files` directory and type this:
```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/bams_mapped_to_XLv10_concatscaf_no_readgroups/* .
```
We also need two more things. The first one is just a list of bam files. You can use the one I made:
```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/angsd/bam.filelist .
```


The second is a phenotype file that tells the software what the phenotype is for each sample: `bin_sex.ybin`; In this case there are two phenotypes: male 0, or female 1.

```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/angsd/bin_sex.ybin .
```

Now we can execute angsd (in a screen) like this:
```
angsd -yBin bin_sex.ybin -doAsso 1 -GL 1 -out out_additive_F1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minInd 4 -bam bam.filelist -P 5 -doCounts 1 -setMinDepthInd 2 -setMaxDepthInd 100 -Pvalue 1
```

This also takes too long to run during class, so let's look at an output file I made earlier:
```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/angsd/pygm_out_additive_F1.lrt0.gz .
zmore pygm_out_additive_F1.lrt0.gz 
```
[pygm_lab_manhattan_angsd.pdf](https://github.com/evansbenj/2024_BIO722/files/14614798/pygm_lab_manhattan_angsd.pdf)
![alt text](https://github.com/evansbenj/2024_BIO722/files/14614798/pygm_lab_manhattan_angsd.pdf))

## Problem 5

With this information in hand, you should be able to make a [Manhattan plot](https://en.wikipedia.org/wiki/Manhattan_plot). Ben will show you some he made.


Now lets quickly look into one way to generate genotypes from our bam files [using bcftools](https://github.com/evansbenj/2024_BIO722/blob/master/5_genotyping_with_bcftools.md)



