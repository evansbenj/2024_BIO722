# Using Stacks with a reference genome

(Or you can go back to using a bash script [here](https://github.com/evansbenj/BIO720/blob/master/3_Lecture_3_Automating_alignment_with_bash.md)).

## Stacks pipeline and setting up a run

[Stacks](http://catchenlab.life.illinois.edu/stacks/manual/) us a software suite that compiles reduced representation genome sequencing data first witih each individual in your analysis, and then across individuals.  It can calculate summary statistics from your data and also output your results in useful formats that other programs can take as input.

We have already used the `process_radtags` module of `Stacks` to de-multiplexing our Illumina data, trimming off linker sequences, and filter sequences based on quality (the proportion of Ns in a read).  There is a very nice online manual [here](http://catchenlab.life.illinois.edu/stacks/manual/) and we will only scratch the surface of what this software can do.

The portion of Stacks that we will use consists of three main steps:
  1. For each individual, sort the data into 'loci' that have one or two alleles (i.e. that are homozygous or heterozygous, respectively).  For analyses with a reference genome, this is accomplished using the program `pstacks` using a `.bam` file as input.  Alternatively if you lack a reference genome, you can use `ustacks` instead of `pstacks`.
  2. Across all individuals, generate a catalog of loci that is a comprehensive list of all genomic regions that have data from at least one individual.  This is accomplished with the `cstacks` program.
  3. Once a catelog of all loci is made, we can compile the data for all individuals to generate a multi-individual genotype for each locus.  This is done with the `sstacks` program

Once loci are compiled within and across individuals, we can use the program `populations` within `Stacks` to analyze the data, including calculating pairwise nucleotide diversity. We can also output the data in different formats that can be analyzed with other software such as [`Phylip`](http://evolution.genetics.washington.edu/phylip/getme.html) and [`Structure`](http://pritchardlab.stanford.edu/structure.html). 

The pipeline of programs within `Stacks` can be run in a batch using a `Perl` script that comes with the program called `refmap.pl`.  This script functions in a similar way to the bash scripts we have used already but it has some added features, such as allowing options to be specified using flags.

To get started, lets first make a directory within the `~/my_monkey_data` directory that has our data called `Stacks_Results`.  To do this, please type this command:

`mkdir ~/my_monkey_data/Stacks_Results`

To make things easier, please also make sure that your bam files from the full data alignment to the rhesus genome are in the `~/my_monkey_data` directory.

## Analysis of population structure

Let's first examine whether population structure is present within our sample by defining two populations.  This can be done in `Stacks` by making a file called a `population_map` in your home directory.  Use your favorite text editor to generate a file that contans this tab delimited information:

```
PF515_sorted	population_1
PM561_sorted	population_1
PM565_sorted	population_1
PM566_sorted	population_1
PM567_sorted	population_1
PM582_sorted	population_2
PM584_sorted	population_2
PM592_sorted	population_2
PM602_sorted	population_2
```

Note that the term tab-delimited means that there is a tab between the columns of information. Please make sure the columns have a tab between them (if you copy and paste the text above it probably will not have a tab). This file will be used to tell Stacks that the first five samples are from one population and the last four samples are from another population. In the example population map above the file names must match the names of each of your bam files, but without the `.bam` suffix.

One way to quantify population structure is using the F-statistic (F<sub>ST</sub>).  F<sub>ST</sub> is an index of population structure that ranges from zero (no population structure) to one (two populations are each fixed for different alleles.  Let's calulate F<sub>ST</sub> between the two populations specified avove using `Stacks`.  To do this, please type:

```
/usr/local/stacks/bin/ref_map.pl -S -b 1 -n 0 \
	-O ~/population_map \
	-o ~/my_monkey_data/Stacks_Results \
   	-s ~/my_monkey_data/PF515_sorted.bam \
    -s ~/my_monkey_data/PM561_sorted.bam \
    -s ~/my_monkey_data/PM565_sorted.bam \
    -s ~/my_monkey_data/PM566_sorted.bam \
    -s ~/my_monkey_data/PM567_sorted.bam \
    -s ~/my_monkey_data/PM582_sorted.bam \
    -s ~/my_monkey_data/PM584_sorted.bam \
    -s ~/my_monkey_data/PM592_sorted.bam \
    -s ~/my_monkey_data/PM602_sorted.bam \
   	-e /usr/local/stacks/bin/ -X "populations:--fstat"
```

In this command, the backslashes `\` just indicate that the command is continued on the next line.  The program we are executing is a Perl script caled `ref_map.pl`.  Similar to the bash scripts we wrote earlier, this program just executes a bunch of other prorgams that come in the `Stacks` package. As you can see in the [`Stacks` manual](http://catchenlab.life.illinois.edu/stacks/comp/ref_map.php) the flags `-S`, `-b`, and `-n` respectively tell Stacks to disable recording the data in an SQL database (this is beyond the scope of this workshop), process the batchID 1, and allow zero mismatches between loci when building the catalog.  Other flags such as `-O`, `-o`, and `-s` respectively tell `Stacks` where the population map file is, where to write the results, and where the input bam files for each individual are.  The `-e` flag tells `Stacks` where the binary files that `refmap.pl` executes are located.

We can also pass some of the programs referenced by `refmap.pl` some additional commands using the `-X` flag.  Here we have used this flag at the end to pass the program `populations` a this flag: `--fstat`, which tells the program `populations` to calculate F<sub>ST</sub> using the population map that we specified using the `-O` flag.  We have additionally specified a directory to write our results to using the `-o` flag.  The `-e` flag tells the computer where the executable files that are referenced by `ref_map.pl` are (these are programs such as `cstacks` and `populations`).

If you now go to the `Stacks_Results` directory (`cd ~/my_monkey_data/Stacks_Results`) and list the files in this directory (`ls`) you should see a bunch of compressed files that have a `gz` suffix.  For each sample we have a file whose name includes the word `alleles`, one with `matches`, one with `snps`, and one with `tags`.  Some details of the contents of these files is available in on [this page](http://catchenlab.life.illinois.edu/stacks/manual/index.php#files) of the Stacks manual.  The file called `batch_1.fst_population_1-population_2.tsv` has the results of the F<sub>ST</sub> calculation for each variable position. This was calculared as we requested by the program `populations` within the `Stacks` software suite.  We can look at this by typing:

`more batch_1.fst_population_1-population_2.tsv`

A summary across all sites is available in another file called `batch_1.fst_summary.tsv`.  Please check this out like this:

`more batch_1.fst_summary.tsv`

## Problem 5

With this information in hand, you should be able to calculate F<sub>ST</sub> for a different comparison by modifying your population_map file.  Please calculate F<sub>ST</sub> between one population comprised of samples PM565, PM566, and PM 567 and the other population comprised of the other samples.

## Now let's move on to [learn more about Stacks](https://github.com/evansbenj/BIO720/blob/master/5_more_on_Stacks.md).



