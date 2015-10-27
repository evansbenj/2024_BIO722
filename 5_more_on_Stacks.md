# Calculating pairwise nucleotide diversity with Stacks

(Or you can go back to [Using Stacks with a Reference Genome](https://github.com/evansbenj/BIO720/blob/master/4_Using_Stacks_to_analyze_your_bam_files.md)).

If we want to calculate pairwise nucleotide diversity on the chromosome to which we mapped our data, we can use the `populations` module of `Stacks` as follows:

`/usr/local/stacks/bin/populations -P path_to_Stacks_Results_folder/Stacks_Results -b 1 -r 1.0 -t 36`

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

Ben will now collect the data you collected, summarize it, compare it to the diversity and divergence estimates from the X, and discuss what implications it has for social system evolution in our monkeys.

# Using a whitelists and blacklists to explore polymorphism inside and outside of genes

One concern about these statistics is that genomic regions that contain genes may be subject to effects of natural selection. To test this, we can use annotation information for our reference genome to identify and exclude RADtags that are located within genes. To do this we will download annotations for the rhemac2 genome and use a perl script to generate a "blacklist" of genomic regions to exclude from our analysis.

## Downloading an annotation of rhesus genes

One advantage to working with an annotated genome sequence is that we can easily access information about gene locations. For the rhesus genome we are working with, we can download the coordinates of known genes [here](http://genome.ucsc.edu/cgi-bin/hgTables). Please visit this page now.

To download the annotations for rhemac2, please select "Mammal" from the clade menue, "Rhesus" from the genome menu, and "Jan. 2006 (MGSC Merged 1.0/rhemac2)" from the assembly menu.  In the group menu please select "Genes and Gene predictions" and for the track menu please select "RefSeq Genes".  In the table menu, please select "refGene" and in the region option please click the "genome" button.  Now change the output format to "BED - browser extensible data" and enter in the output file a name for your annotation file such as "rhemac2.bed". Once you have done all this, please click the "get output" button. Please upload this to your info account.

You can check out the contents of this file by typing:

`more rhemac2.bed`

This [page](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) provides information about many genomic formats, including `bed` files. Please check it out.  As you can see, the first three columns are the ones we care about: these tell us the chromosome number, start and end coordinates of genes.  Please note that the numbering format for bed files is somewhat weird; the start coordinate is based on a numbering system beginning with zero, not one. And the end coordinate refers to the base after the end of the feature. As stated on this page, the first 100 base pairs of a chromosome would have bed coordinates of 0 to 100, and would span bases that are numbered 0 to 99.

## Making a whitelist from the bed file for use with `Stacks`

I (hopefully) wrote a perl script that hopefully will generate a list of RAD loci from the `batch_1.catalog.tags.tsv` file that fall within the intervals provided in the `rhemac2.bed` file. We can then use this to generate polymorphism statistics inside and outside of genes with the `populations` module of `Stacks`.




# Now let's use `Stacks` to make a phylogenetic tree [here](https://github.com/evansbenj/BIO720/blob/master/6_Making_a_phylogenetic_tree_with_Stacks.md).
