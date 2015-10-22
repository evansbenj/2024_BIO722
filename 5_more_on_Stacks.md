# More on Stacks:  Whitelists, blacklists, using individual modules, and summary statistics

(Or you can go back to [Using Stacks with a Reference Genome](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/6_Using_Stacks_with_a_reference_genome.md)).

## Whitelists, Blacklists, and Calculating Summary Statistics with `Populations`

Now lets try to calculate a summary statistic from specific genomic regions using the `populations` module of stacks. To accomplish this, lets work with a larger dataset that I made earlier.  Please switch to this directory:

`cd /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results`

Let's first look at the catalog file like this:

`more batch_1.catalog.tags.tsv`

The first line begins with a hash (`#`) symbol and is reserved for comments.  The next lines have columns of text.  The 4th column lists the chromosome number and position of each tag.  We are going to generate a file in which we sample random SNPs from 1000 tags but we want to exclude data from the X chromosome because there are differences in copy number between males and females (i.e. two in XX females and one in XY males).  In order to do this, we can create a list of tags to include (a `whitelist`) or exclude (a `blacklist`), which is just a list of the numbers in the 3rd column that correspond with the `chrX` in the 4th column.  To generate a list based on this criterion, please use this `unix` command:

`awk '$4 ~ /chrX/ {print $3}' batch_1.catalog.tags.tsv > ~/chrX_list`

This uses a `unix` function called `awk`.  It says to print the number in column 3 to a file in your home directory (`~`) called `chrX_list` whenever the value in column 4 is equal to `chrX`.

We can view this file by typing:

`more ~/chrX_list`

Now, if we want to calculate pairwise nucleotide diversity on the autosomes, we can use the `populations` module of `Stacks` with the `blacklist` flag as follows:

`/apps/stacks/1.29/bin/populations -P ~/monkey/Stacks_Results -b 1 -r 1 -t 36 -B ~/chrX_list`

This command directs the `populations` module of `Stacks` to write results to a directory specified by the `-P` flag. The `-r` flag says we want to only print data where 100% of the individuals have a genotype.  The `-b` and `-t` flags specify, respectively, that `populations` should focus on batch_ID number 1 (you can work with special IDs if you need to but this is beyond the scope of this workshop) and that `populations` should use 36 threads to do the calculations.  The calculations will include all positions except those on the X chromosome, which is blacklisted by the `-B` flag.

We can view the pairwise nucleotide diversity statistic (Pi) in a file called `batch_1.sumstats_summary.tsv` like this:

`more ~/monkey/Stacks_Results/batch_1.sumstats_summary.tsv`

Similarly, if we want to calculate Pi on only the X chromosome, we can use the `chrX_list` as a `whitelist` as follows:

`/apps/stacks/1.29/bin/populations -P ~/monkey/Stacks_Results -b 1 -r 1 -t 36 -W ~/chrX_list`

There are problems with this calculation using only the X chromosome because some of the individuals are males; Ben will discuss these problems and a workaround.

# Now let's use `Stacks` to make a `Structure` input file [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/8_Stacks_and_Structure.md).
