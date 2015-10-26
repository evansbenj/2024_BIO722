# Calculating pairwise nucleotide diversity with Stacks

(Or you can go back to [Using Stacks with a Reference Genome](https://github.com/evansbenj/BIO720/blob/master/4_Using_Stacks_to_analyze_your_bam_files.md)).

If we want to calculate pairwise nucleotide diversity on the chromosome to which we mapped our data, we can use the `populations` module of `Stacks` as follows:

`/usr/local/stacks/bin/populations -P path_to_Stacks_Results_folder/Stacks_Results -b 1 -r 100 -t 36`

As detailed in the [manual](http://catchenlab.life.illinois.edu/stacks/comp/populations.php), this command directs the `populations` module of `Stacks` to write results to a directory specified by the `-P` flag. The `-r` flag says we want to only print data where 100% of the individuals have a genotype.  The `-b` and `-t` flags specify, respectively, that `populations` should focus on batch_ID number 1 (you can work with special IDs if you need to but this is beyond the scope of this workshop) and that `populations` should use 36 threads to do the calculations.  

We can view the pairwise nucleotide diversity statistic (Pi) in a file called `batch_1.sumstats_summary.tsv` like this:

`more ~/monkey/Stacks_Results/batch_1.sumstats_summary.tsv`

As indicated in the header, the pairwise nucleotide diversity statistic (Pi) should be listed in the sixth to last column on the bottom row. For autosomes this statistic is accurate but there are problems with this calculation using only the X chromosome because some of the individuals are males; Ben will discuss these problems and a workaround he did with a perl script. The estimate of pi he recovered from the Stacks analysis of the bam files is 0.00108 substutitons per site.

Now lets take some time for you to all analyze your bam files. Please let Ben know what your calculation of pi is for each of your chromosomes.

# Calculating divergence to the reference genome using `POPBAM`

As some of you may realize, the mutation rates in males and females are not equivalent. In many species, including all primates, more cell divisions occur to make sperm than eggs, and for this reason (we think) the mutation rate is higher in males than in females.  Because X chromosomes spend 2/3rds of their time in females, we expect lower polymorphism on the X by virtue of this.  We can control for differences in mutation rate very simply – by taking the ratio of pi to divergence in the X and the aDNA. We calculated pi above, now we need divergence.

[`POPBAM`](http://popbam.sourceforge.net/) is a program that is useful for calculating several population genetic statistics. Some of the functionality overlaps with that of `Stacks`.  We will use `POPBAM` to calculate the level of divergence of one of our samples to the reference genome.

First we need to add some additional information to the readgroup header of our `bam` file.  Please type this:

```
samtools view -H XXX.bam > header.sam
perl -pi.old -e 's{PL:illumina}{PL:illumina\tPO:POP1}g' header.sam
samtools reheader header.sam XXX.bam > XXX_new.bam
samtools index XXX_new.bam
```

The first line uses samtools to make a text file called `header.sam` that contains the header information for the file called `XXX.bam`.  You should use one of your sorted bam files for this. The second line uses `Perl` to search and replace text in the header.sam file.  Basically this adds text to the readgroup portion of this file. Thie third line uses samtools to change the header of our file. And the fourth line makes a new index file for our new bam file.

Now we can get divergence information using `POPBAM` like this:

`popbam diverge -o 0 XXX_new.bam chrX -f path_to_reference_chromosome/chrXXX.fa chrXXX > divergence.txt`

Don't forget to specify the chromosome that you are surveying at the end (`chrXXX`). This will write the output of `POPBAM` to a text file called `divergence.txt`.  The divergence Ben got for PF515 was 0.00524.

# Now let's use `Stacks` to make a phylogenetic tree [here](https://github.com/evansbenj/BIO720/blob/master/6_Making_a_phylogenetic_tree_with_Stacks.md).
