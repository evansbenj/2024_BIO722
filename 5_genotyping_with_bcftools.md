# Genotyping with Bcftools

Or you can go back to [using angsd for GWAS](4_Using_angsd_for_GWAS.md)

We just made some symbolic links to lots of bam for angsd analysis. We can also use these files as input for genotyping.

A common way that genotype information is conveyed is the `variant call format` – `.vcf`, which is is described [here](https://en.wikipedia.org/wiki/Variant_Call_Format). Another new format introduced by the [`Genome Analysis Toolkit`](https://en.wikipedia.org/wiki/Variant_Call_Format) of the Broad Institute is called the genomic variant call format – `.gvcf`; this is described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format).

Below is a command that you can use to make a `.vcf` file for one chromosome for one sample. For fun (if you and some free time) you can write a bash script to make a `.vcf` file for all chromosomes and for each `.bam` file using a nested loop.

Don't run this now because it will take too long to finish in class.  But here is an example
```
bcftools mpileup -f REFERENCE sample1.bam sample2.bam .. | bcftools call -V indels --format-fields GQ -m -O z -o OUTPUTFILE.vcf.gz
```
For example:
```
bcftools mpileup -f ../reference/XENLA_10.1_genome.fa.gz pyg_mal_Z23368_TACAT_sorted.bam pyg_mal_Z23376_CTTCCA_sorted.bam | bcftools call -V indels --format-fields GQ -m -O z -o allsamples_merged_sorted.bam.vcf.gz
```

If you wanted to run this, you should do it using `screen` as we discussed earlier.

The control above first uses the samtools mpileup command. You can check out what the options are by typing `./samtools mpileup`.  You will see that the -d command asks `samtools` to allow very deep depth of covergae for each position (this is probably overkill). The `-ugf` commands tell `samtools` respectively to generate uncompressed VCF/BCF output, generate genotype likelihoods in BCF format, and that the reference file is in fasta format and that it has a faidx index. The `-t` flag tells `samtools` what information to output for each genotype (total depth and per allele depth in this case). Then the sorted bam files to genotype are listed. `samtools` will compute the likelihood of the data given each possible genotype for each position for each sample.  It does not call the variants though. 

This information is then piped to a program called `bcftools` which applies a prior probability to each genotype for each position across all individuals and does the genotype calling. The flags tell `bcftools` to skip indels, output a compressed file, and include genotype qualities in the output. 

Eventually this will finish. In the meantime, we can check out some files that I made earlier. Please make a directory called `vcf_files` and enter this directory. Now please make a symbolic link to this file like this:

```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/vcfs_mapped_to_XLv10_concatscaf/* .
```

And now check out this file like this:

```
more pyg_mal_Z23376_CTTCCA_sorted.bam_rg.bam_Chr9_10S.g.vcf

```

If the vcf files were compressed (.g.vcf.gz) thatn you would need to use `zmore` instead of `more`.

This file was actually made with GATK instead of samtools, but the overall format is similar (GATK gives more information and you can ask bcftools to do this if you want). Depending on time, Ben will discuss other genotyping approaches using GATK, base recalibration, genotype filtering steps, and also some other useful tasks you can perform using vcf files - such as principal components analysis.

OK now let's discuss some [tips for using the Digital Research Alliance of Canada](https://github.com/evansbenj/2024_BIO722/blob/master/6_tips_for_using_the_digital_research_alliance_of_canada.md)
