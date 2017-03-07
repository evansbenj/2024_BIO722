# Genotyping with Samtools and bcftools

(Or you can go back to using `Stacks` to make a Structure files [here](https://github.com/evansbenj/BIO720/blob/master/7_Stacks_and_Structure.md)).

Now lets make a genotype (vcf) file from our bam files.

First, from within this directory: `/2/scratch/XXXXX/full_data_mapped_to_chr9`, please index all of the chr9 files with this script (don't forget to update the `XXXXX` in the path with your username):

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

