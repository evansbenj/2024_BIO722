# Genotyping with Samtools and bcftools

(Or you can go back to using `angsd` for GWAS analysis [here]([https://github.com/evansbenj/BIO720/blob/master/7_Stacks_and_Structure.md](https://github.com/evansbenj/2024_BIO722/blob/master/4_Using_angsd_for_GWAS.md)).

Remember previously we made some symbolic links to some sorted bam files I made to chr9 in a directory `/2/scratch/ZZZ/full_data_mapped_to_chr9`.  I made these so we could work with them in class (they take a while to generate). Please enter this directory, and check that they are there.  If not, please make them again like this:
```
ln -s /1/scratch/ben/PF515_chr9_sorted.bam
ln -s /1/scratch/ben/PM561_chr9_sorted.bam
ln -s /1/scratch/ben/PM565_chr9_sorted.bam
ln -s /1/scratch/ben/PM566_chr9_sorted.bam
ln -s /1/scratch/ben/PM567_chr9_sorted.bam
ln -s /1/scratch/ben/PM582_chr9_sorted.bam
ln -s /1/scratch/ben/PM584_chr9_sorted.bam
ln -s /1/scratch/ben/PM592_chr9_sorted.bam
ln -s /1/scratch/ben/PM602_chr9_sorted.bam
```

A common way that genotype information is conveyed is the `variant call format` – vcf, which is is described [here](https://en.wikipedia.org/wiki/Variant_Call_Format). Another new format introduced by the [`Genome Analysis Toolkit`](https://software.broadinstitute.org/gatk/) of the Broad Institute is called the genomic variant call format – gvcf; this is described [here](http://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf).

Now lets make a vcf file from our bam files.

First, please index all of the chr9 files with this script. You will need to copy and paste this script into a text editor, save it under some name, and don't forget to update the `XXXXX` in the `path_to_data` variable with your username.  Then you can execute it.

```
#!/bin/bash                                                                                            

path_to_data="/2/scratch/XXXXX/full_data_mapped_to_chr9"
chromosome="chr9"

individuals="PF515                                                                                     
PM561                                                                                                  
PM565                                                                                                  
PM566                                                                                                  
PM567                                                                                                  
PM582                                                                                                  
PM584                                                                                                  
PM592                                                                                                  
PM602"

for each_individual in $individuals
do
    echo ${each_individual}
    samtools index $path_to_data/${each_individual}_${chromosome}_sorted.bam
done
```

Now let's make some symbolic links to some recently updated software. The reason we are doing is is that the newest release of samtools has some new functionality that we will need:

```
ln -s /home/ben/samtools_2016/bin/samtools
ln -s /home/ben/samtools_2016/bcftools-1.3.1/bcftools
```


Now please make a directory called `my_chr9`, enter this directory, and make some more symbolic links to some chr9 reference files. We previously made these for a different chromosome using commands we went through earlier and I used the same commands to make these files for chromosome 9.  We're making links to these files now to save some time:

```
ln -s /1/scratch/ben/chr9_reference_genome/chr9.dict
ln -s /1/scratch/ben/chr9_reference_genome/chr9.fa
ln -s /1/scratch/ben/chr9_reference_genome/chr9.fa.amb
ln -s /1/scratch/ben/chr9_reference_genome/chr9.fa.ann
ln -s /1/scratch/ben/chr9_reference_genome/chr9.fa.bwt
ln -s /1/scratch/ben/chr9_reference_genome/chr9.fa.fai
ln -s /1/scratch/ben/chr9_reference_genome/chr9.fa.pac
ln -s /1/scratch/ben/chr9_reference_genome/chr9.fa.sa
```


Now, we have everything we need to make a genotype (vcf) file with all of the samples, including:
* sorted bam files for each sample
* an index for each bam file (.bai file)
* a genome (actually a chromosome in our case - chr9.fa) and also supporting files (chr9.fa.dict, chr9.fa.fai, etc).

Please start a screen and type this command:

```
./samtools mpileup -d8000 -ugf ./my_chr9/chr9.fa -t DP,AD PF515_chr9_sorted.bam PM561_chr9_sorted.bam PM565_chr9_sorted.bam PM566_chr9_sorted.bam PM567_chr9_sorted.bam PM582_chr9_sorted.bam PM584_chr9_sorted.bam PM592_chr9_sorted.bam PM602_chr9_sorted.bam | ./bcftools call -V indels --format-fields GQ -m -O z -o allsamples_chr9_merged_sorted.bam.vcf.gz

```
This will take a while so please exit the screen by typing `ctrl-a` and `ctrl-d`

The control above first uses the samtools mpileup command. You can check out what the options are by typing `./samtools mpileup`.  You will see that the -d command asks `samtools` to allow very deep depth of covergae for each position (this is probably overkill). The `-ugf` commands tell `samtools` respectively to generate uncompressed VCF/BCF output, generate genotype likelihoods in BCF format, and that the reference file is in fasta format and that it has a faidx index. The `-t` flag tells `samtools` what information to output for each genotype (total depth and per allele depth in this case). Then the sorted bam files to genotype are listed. `samtools` will compute the likelihood of the data given each possible genotype for each position for each sample.  It does not call the variants though. 

This information is then piped to a program called `bcftools` which applies a prior probability to each genotype for each position across all individuals and does the genotype calling. The flags tell `bcftools` to skip indels, output a compressed file, and include genotype qualities in the output. 

Eventually this will finish. In the meantime, we can check out a file I made earlier.  Please make a symbolic link to this file like this:

```
ln -s /1/scratch/ben/allsamples_chr9_merged_sorted.bam.vcf.gz
```

And now check out this file like this:

```
zmore allsamples_chr9_merged_sorted.bam.vcf.gz

```

Here you need to use `zmore` instead of `more` because the file is compressed (gzipped).

# Making a tab delimited genotype file

The information in the vcf file is rich but also sometimes TMI. I usually convert vcf files to tab delimited files that I use as input for millions of perl scripts.  You can co this using a program called `vcftools` like this:

```
/usr/local/tabix/tabix -f -p vcf allsamples_chr9_merged_sorted.bam.vcf.gz
zcat allsamples_chr9_merged_sorted.bam.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > allsamples_chr9_merged_sorted.bam.vcf.gz.tab
```
 
 This takes a very little while. You can check out the tab file you made like this:
 
 ```
 more allsamples_chr9_merged_sorted.bam.vcf.gz.tab
 ```

Depending on time, Ben will discuss other genotyping approaches using GATK, base recalibration, genotype filtering steps, and also some other useful tasks you can perform using vcf files - such as principal components analysis.
