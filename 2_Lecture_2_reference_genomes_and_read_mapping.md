# BJE Lecture 2

(or you can go back to the de-multipexing page [here](https://github.com/evansbenj/BIO720/blob/master/1_Lecture_1.md)).

Depending on your organism of study, there may or may not be a relatively closely related genome sequence to work with.  Depending on your research question, this may or may not be useful.  In our example study on Tonkean macaques, we are  interested in quantifying molecular polymorphism on the X chromosome and comparing it to polymorphism on the autosomes.  For this reason, the genomic location of the data is crucial and we can benefit from the complete genome sequence of a closely related species of macaque monkey, the rhesus macaque (*Macaca mulatta*).  We will use a program called [`bwa`] (http://sourceforge.net/projects/bio-bwa/files) and also [`samtools`](http://samtools.sourceforge.net/), to map our data to individual chromosomes of the rhesus macaque.  Normally one would map reads to an entire genome because the data were generated from a complete genome, but in our case we are doing only an example analysis and we will each work on an individual chromosome. Ben will assign each of you a chromosome to work on.

## A note about "completely" sequenced genomes

FYI, essentially all completely sequenced genomes are not in fact completely sequenced.  
- Regions such as centromeric and telomeric regions and some portions of sex-specific sex chromosomes contain many repetitive elements that pose challenges to sequencing and assembly.  
- Sometimes the individual sequenced is female, so no Y chromosome is available.  Sometimes (usually?) when a genome is said to be "complete" it actually is a bunch of "contigs", or contiguous sequence, that may or may not be assembled into "scaffolds" that include contigs plus Ns to represent intervenining regions that are not yet sequenced.  
- And even then, we can expect sequence and assembly errors in our reference genome that make it different from the real genome sequence.  
- On top of that, there is population level variation to contend with, including SNPs and insertion deletion events.  This makes our samples different from any reference genome as well.

As an example, let's look at some information on the "completely" sequenced genomes of [some frogs](http://www.xenbase.org/other/static/ftpDatafiles.jsp).  Of interest is the N50 statistic of a genome assembly, which is defined [here](https://en.wikipedia.org/wiki/N50_statistic).

## Preparing your reference genome

Reference genomes for many sequences are available at multiple publicly available databases.  We can download the complete genome sequence for the rhesus macaque from the [USC genome browser](http://hgdownload.cse.ucsc.edu/downloads.html#rhesus).  I did this earlier because it takes a while.  The whole genome comes as a fasta-formatted file, and I split it up into individual fasta files corresponding with each of the chromosomes.  These are located in this directory:

`/1/scratch/BIO720_Bens_section/rhesus_chromosomes/`

Now check out what is in this directory by typing this:

`ls /1/scratch/BIO720_Bens_section/rhesus_chromosomes/`

Ben will assign you a chromosome to work with.  Please copy this chromosome sequence to a new directory that you make like this:

`mkdir ~/my_monkey_chromosome`

and then copy your chromosome to that directory like this:

`cp /1/scratch/BIO720_Bens_section/rhesus_chromosomes/chrXXX.fa ~/my_monkey_chromosome/.` 

Here and henceforth, you will need to change the `chrXXX.fa` part to match whatever chromosome Ben assigned to you.  For example, if you are working on chromosome 9, you should type this:

`cp /1/scratch/BIO720_Bens_section/rhesus_chromosomes/chr9.fa ~/my_monkey_chromosome/.` 

You can check whether you succeeded and also see how big this datafile is like this:

`ls -l ~/my_monkey_chromosome/chrXXX.fa` 

Before we map our data to this reference genome, we need to generate some files that will be used in the mapping process.  This can be done in three steps:

1. `bwa index -a bwtsw ~/my_monkey_chromosome/chrXXX.fa`

  The `/apps/bwa/0.7.12/bwa` command tells the computer to execute the bwa program.  The `index` command tells `bwa` to generate index files from the rhesus genome file that is indicated by the `/home/datasets/2015_Ben_Evans/rhesus_chromosomes/chrXXX.fa`. The `-a bwtsw` flag specifies the indexing algorithm for `bwa` to use.  
  
  (But delete the two XXs - I added them in to prevent people from just copying and pasting things)
  
  This step will take a few minutes.

2. We now need to to generate another file using `samtools`.  Please type this:

  `samtools faidx ~/my_monkey_chromosome/chrXXX.fa`

  Here, the `samtools` command tells the computer to execute the `samtools` program.  The `faidx` option tells samtools to generate a file called `chrXXX.fai` in which each line has information for one the contigs within the reference genome, including the contig name, size, location and other information.  Our reference genome has a contig for each chromosome.

3.  The third thing we need to do is to generate a `.dict` file with a program called [`picard`](http://broadinstitute.github.io/picard/).  To do this, first enter the `my_monkey_chromosome` directory like this:

`cd ~/my_monkey_chromosome`

Once you are in, please type this command:

  `java -jar /usr/local/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=chrXXX.fa OUTPUT=chrXXX.dict`

  As before, you will need to change the `chrXXX` in this command to match the chromosome you are working with.  This should generate a file called `~/my_monkey_chromosome/chrXXX.dict`

## Mapping the data to the reference genome

Now we can align the data from each individual to the reference genome using [`bwa`] (http://sourceforge.net/projects/bio-bwa/files). Because this takes a while to do with the fill dataset, I made a datasubset file that we can work with in class. You can find this here:  

`/1/scratch/BIO720_Bens_section/subset_data/PF515_subset.fastq`

Please copy this to a new directory that you can make like this:

`mkdir ~/my_monkey_data`

and then copy it to this directory like this:

`cp /1/scratch/BIO720_Bens_section/subset_data/PF515_subset.fastq ~/my_monkey_data`

Now let's map this data subset from one individual to the reference genome using `bwa` as follows:

`bwa aln reference_genome data.fastq > data.sai`

For example, for this individual (PF515) you could type this

`bwa aln ~/my_monkey_chromosome/chrXXX.fa ~/my_monkey_data/PF515_subset.fastq > ~/my_monkey_data/PF515_subset.sai`

(but with the `chrXXX.fa` changed to match the chromosome you are working on.)

As you can see in the [`bwa` manual] (http://bio-bwa.sourceforge.net/bwa.shtml), the `aln` flag tells `bwa` to align the reads to the reference genome (or, using the `bwa` jargon, find the coordinates of the reference genome that map each read. There are several additional options you could specify if you want.  For example, you could use the `-M` flag to set a limit on the number of mismatches between a read and the reference genome.  Because we are mapping data from one species to a reference genome from another, we will not do this.

This command generates an intermediate file with the `.sai` suffix (which stands for `suffix array index`). Now we need to generate a `.sam` formatted file from our `.sai` files.  A `.sam` file is a tab delimited text file that contains the alignment data.  The format for this command is:

`bwa samse reference_genome.fa data.sai data.fastq > data.sam`

We also need to add a header to each `.sam` file, so we can type this command:

`bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:XXX\tPL:illumina" reference.fa data.sai data.fastq > data.sam`

Bute here you need to change the "XXX" to match the sample you are working with (for example change XXX to PF515 if you are working with sample PF515).

Now we can generate a `.bam` file.  A `.bam` formatted file is a binary version of the `.sam` file.

`samtools view -bt reference_genome -o data.bam data.sam`

Sort the `.bam` file:

`samtools sort data.bam data_sorted`

Make an index for the bam file, which is a `.bai` file:

`samtools index data_sorted.bam`

## Problem 3: Assessing coverage

Samtools can provide information on the number of reads for each position of the reference sequence for which there are data.  You can see this information by typing this:

`samtools depth XXX_sorted.bam`

Where `XXX` is the sample ID number.  If you want to know the average depth across all sites, you could type this:

`samtools depth XXX_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'`

Here the vertical bar `|` is a "pipe" that sends the information from the command before it to the command after it.  So the data you generated from `samtools` will be parsed with the unix `awk` command.  This will add the values of the third column `$3` to a variable called `sum` and then at the end (`END`) print out the word `Average` followed by the quotient `sum/NR` where `NR` is a built in variable that keeps track of the number of records.  A good description of `awk` is [here](http://www.folkstalk.com/2011/12/good-examples-of-awk-command-in-unix.html).

Now it is your turn.  Using the manual for [samtools](http://www.htslib.org/doc/samtools-0.1.19.html) please quantify how many reads mapped to your chromosome for each individual and how many failed to map.  Do you know why so many failed to map?

## Practice Problem

Please write a bash script that will provide the coverage statistic for all bam files in a folder that have the suffix "_sorted.bam".  You can (and should) use the internet for tips.

## OK, if this all went smoothly we are now ready to automate the alignments with a bash script.  Please click [here](https://github.com/evansbenj/BIO720/blob/master/3_Lecture_3_Automating_alignment_with_bash.md) to go to the next page.
