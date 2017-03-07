# Genotyping with Samtools and bcftools

(Or you can go back to using `Stacks` to make a Structure files [here](https://github.com/evansbenj/BIO720/blob/master/7_Stacks_and_Structure.md)).

Now lets make a genotype (vcf) file from our bam files.

First, please index all of the chr9 files with this script. You will need to copy and paste this script into a text editor, save it under some name, and don't forget to update the `XXXXX` in the path with your username.  Then you can execute it.

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

Let's make some symbolic links to some recently updated software:

```
ln -s /home/ben/samtools_2016/bin/samtools
ln -s /home/ben/samtools_2016/bcftools-1.3.1/bcftools
```

Now, using these recently updated versions, we can make a genotype (vcf) file with all of the samples like this:

```
./samtools mpileup -d8000 -ugf REFERENCE.fa -t DP,AD PF515_chr9_sorted.bam PM561_chr9_sorted.bam PM565_chr9_sorted.bam PM566_chr9_sorted.bam PM567_chr9_sorted.bam PM582_chr9_sorted.bam PM584_chr9_sorted.bam PM592_chr9_sorted.bam PM602_chr9_sorted.bam | ./bcftools call -V indels --format-fields GQ -m -O z -O z -o allsamples_chr9_merged_sorted.bam.vcf.gz

```

