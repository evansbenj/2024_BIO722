# Calculating pairwise nucleotide diversity with Stacks

(Or you can go back to [Using Stacks with a Reference Genome](https://github.com/evansbenj/BIO720/blob/master/4_Using_Stacks_to_analyze_your_bam_files.md)).

If we want to calculate pairwise nucleotide diversity on the chromosome to which we mapped our data, we can use the `populations` module of `Stacks` as follows:

`/usr/local/stacks/bin/populations -P ~/monkey/Stacks_Results -b 1 -r 1 -t 36`

As detailed in the [manual](http://catchenlab.life.illinois.edu/stacks/comp/populations.php), this command directs the `populations` module of `Stacks` to write results to a directory specified by the `-P` flag. The `-r` flag says we want to only print data where 100% of the individuals have a genotype.  The `-b` and `-t` flags specify, respectively, that `populations` should focus on batch_ID number 1 (you can work with special IDs if you need to but this is beyond the scope of this workshop) and that `populations` should use 36 threads to do the calculations.  

We can view the pairwise nucleotide diversity statistic (Pi) in a file called `batch_1.sumstats_summary.tsv` like this:

`more ~/monkey/Stacks_Results/batch_1.sumstats_summary.tsv`

As indicated in the header, the pairwise nucleotide diversity statistic (Pi) should be listed in the sixth to last column on the bottom row. For autosomes this statistic is accurate but there are problems with this calculation using only the X chromosome because some of the individuals are males; Ben will discuss these problems and a workaround.

# Calculating divergence to the reference genome using `POPBAM`

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

Don't forget to specify the chromosome that you are surveying at the end (`chrXXX`). This will write the output of `POPBAM` to a text file called `divergence.txt`


# More on Stacks:  Whitelists, blacklists, using individual modules, and summary statistics in `Stacks`

## Whitelists, Blacklists, and Calculating Summary Statistics with `Populations` 

Now lets try to calculate a summary statistic from specific genomic regions using the `populations` module of stacks. To accomplish this, lets work with a larger dataset that I made earlier.  Please copy this directory to your `my_monkey_data` directory like this:

`cp -r /1/scratch/BIO720_Bens_section/Stacks_Results_chrX ~/my_monkey_data/.`

Let's unzip the catalog file and check it out like this:

`gunzip batch_1.catalog.tags.tsv.gz`

and then

`more batch_1.catalog.tags.tsv`

The file should contain multiple columns of text.  The 4th and 5th columns list the chromosome number and chromosome position of each tag respectively.  We are going to generate a file in which we sample random SNPs from 1000 tags but we want to exclude data from the X chromosome because there are differences in copy number between males and females (i.e. two in XX females and one in XY males).  In order to do this, we can create a list of tags to include (a `whitelist`) or exclude (a `blacklist`), which is just a list of the numbers in the 3rd column that correspond with the `chrX` in the 4th column.  To generate a list based on this criterion, please use this `unix` command:

`awk '$4 ~ /chrX/ {print $3}' batch_1.catalog.tags.tsv > ~/chrX_list`

This uses a `unix` function called `awk`.  It says to print the number in column 3 to a file in your home directory (`~`) called `chrX_list` whenever the value in column 4 is equal to `chrX`.

We can view this file by typing:

`more ~/chrX_list`

# Now let's use `Stacks` to make a `Structure` input file [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/8_Stacks_and_Structure.md).
