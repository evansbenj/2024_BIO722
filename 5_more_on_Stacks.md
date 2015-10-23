# Calculating pairwise nucleotide diversity$\pi$ with Stacks

(Or you can go back to [Using Stacks with a Reference Genome](https://github.com/evansbenj/BIO720/blob/master/4_Using_Stacks_to_analyze_your_bam_files.md)).

If we want to calculate pairwise nucleotide diversity on the chromosome to which we mapped our data, we can use the `populations` module of `Stacks` with the `blacklist` flag as follows:

`/usr/local/stacks/bin/populations -P ~/monkey/Stacks_Results -b 1 -r 1 -t 36`

As detailed in the [manual](http://catchenlab.life.illinois.edu/stacks/comp/populations.php), this command directs the `populations` module of `Stacks` to write results to a directory specified by the `-P` flag. The `-r` flag says we want to only print data where 100% of the individuals have a genotype.  The `-b` and `-t` flags specify, respectively, that `populations` should focus on batch_ID number 1 (you can work with special IDs if you need to but this is beyond the scope of this workshop) and that `populations` should use 36 threads to do the calculations.  

We can view the pairwise nucleotide diversity statistic (Pi) in a file called `batch_1.sumstats_summary.tsv` like this:

`more ~/monkey/Stacks_Results/batch_1.sumstats_summary.tsv`

As indicated in the header, the pairwise nucleotide diversity statistic (Pi) should be listed in the sixth to last column on the bottom row. For autosomes this statistic is accurate but there are problems with this calculation using only the X chromosome because some of the individuals are males; Ben will discuss these problems and a workaround.

# More on Stacks:  Whitelists, blacklists, using individual modules, and summary statistics

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
