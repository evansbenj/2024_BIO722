# Aligning to a Reference Genome

(or you can go back to the Quality Control and Trimming page [here](https://github.com/evansbenj/BIO720/blob/master/1_Lecture_1.md)).

Depending on your organism of study, there may or may not be a relatively closely related genome sequence to work with.  Depending on your research question, this may or may not be useful.  In our example study on the frog *Xenopus pygmaeus*, we are  interested in ultimately doing a genome-wide association study to identify the sex determining region.  For this reason, the genomic location of the data is crucial and we can benefit from the complete genome sequence of a closely related species of African clawed frog (*Xenopus laevis*).  We will use a program called [`bwa`](https://bio-bwa.sourceforge.net/) and also [`samtools`](http://samtools.sourceforge.net/), to map our data to individual chromosomes of *Xenopus laevis*. 

## A note about "completely" sequenced genomes

FYI, essentially all completely sequenced genomes are not in fact completely sequenced.  
- Regions such as centromeric and telomeric regions and some portions of sex-specific sex chromosomes contain many repetitive elements that pose challenges to sequencing and assembly.  
- Sometimes the individual sequenced is female, so no Y chromosome is available.  Sometimes (usually?) when a genome is said to be "complete" it actually is a bunch of "contigs", or contiguous sequence, that may or may not be assembled into "scaffolds" that include contigs plus Ns to represent intervenining regions that are not yet sequenced.  
- And even then, we can expect sequence and assembly errors in our reference genome that make it different from the real genome sequence.  
- On top of that, there is population level variation to contend with, including SNPs and insertion/deletion events.  This makes our samples different from any reference genome as well.

As an example, let's look at some information on the "completely" sequenced genomes of [some frogs](http://www.xenbase.org/other/static/ftpDatafiles.jsp).  Of interest is the N50 statistic of a genome assembly, which is defined [here](https://en.wikipedia.org/wiki/N50_statistic).

## Preparing your reference genome

Reference genomes for many sequences are available at multiple publicly available databases.  We can download the complete genome sequence for *Xenopus laevis* from the [Xenbase](https://www.xenbase.org/xenbase/static-xenbase/ftpDatafiles.jsp). 

Please exit out of your fq directory and make a new directory called `Reference_Genome`. Enter that directory. You can download the *Xenopus laevis* genome assembly using `wget` as follows (please don't do this though):
```
wget https://download.xenbase.org/xenbase/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
```
The whole genome comes as a fasta-formatted file. I did this earlier because it takes a while.  Instead of downloading it, please just make a symbolic link to the file I have:

`ln -s /home/ben/2024_BIO722/2022_pygmaeus/2021_XL_v10_refgenome/XENLA_10.1_genome.fa.gz .`

Now check out the first few lines:

`zcat XENLA_10.1_genome.fa.gz | more`


Before we map our data to this reference genome, we need to generate some files that will be used in the mapping process.  This can be done in three steps:

1. Make an index file.   

    The `bwa` command tells the computer to execute the bwa program.  The `index` command tells `bwa` to generate index files from the *X. laevis* genome file:
```
bwa index XENLA_10.1_genome.fa.gz
```  
This step will take a while, and we could do it in the background using `screen`. Ben will do a demo of how to use screen (just watch, no need to do this):

`screen -S make_an_index_file`
  
then type this:
  
`bwa index XENLA_10.1_genome.fa.gz`
  
Then exit the screen by typing `ctrl-a` then `ctrl-d`
  
You can list the screens you have like this:
  
`screen -ls`

and you can return to the screen you started like this:
  
`screen -r make_an_index_file`
  
when it is done, you can exit and then kill the screen like this:
  
`ctrl-a` then `ctrl-d` and then
  
screen -X -S make_an_index_file kill


2. For some downstram applications, we need to to generate another file using `samtools`.  Please type this:

  `samtools faidx XENLA_10.1_genome.fa.gz`

  Here, the `samtools` command tells the computer to execute the `samtools` program.  The `faidx` option tells samtools to generate a file called `XENLA_10.1_genome.fa.gz.fai` in which each line has information for one the contigs within the reference genome, including the contig name, size, location and other information.  Our reference genome has a contig for each chromosome.

3.  The third thing we need to do is to generate a `.dict` file with a program called [`picard`](http://broadinstitute.github.io/picard/). You can run picard like this:

  `java -jar /usr/local-centos6/picard-tools/picard.jar CreateSequenceDictionary REFERENCE=XENLA_10.1_genome.fa.gz OUTPUT=XENLA_10.1_genome.fa.gz.dict`

This should generate a file called `XENLA_10.1_genome.fa.gz.dict`

Because these steps all take time, let's just make symbolic links to the files I already generated:
```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/2021_XL_v10_refgenome/XENLA_10.1_genome.fa.gz.* .
```

## Mapping the data to the reference genome

Please make a directory called `bam_files` and enter this directory. 

Now we can align the data from each individual to the reference genome using [`bwa`](https://bio-bwa.sourceforge.net/). 

```
bwa mem reference forward.R1.fq.gz reverse.R2.fq.gz | samtools view -Shu - | samtools sort - -
o sampleID_sorted.bam`
samtools index sampleID_sorted.bam
```

For example, you could map our subsetted files for individual Z23337 to the *X. laevis* genome like this:

```
bwa mem ../reference/XENLA_10.1_genome.fa.gz ../fq/Z23337_CTCG_R1_subset.fq ../fq/Z23337_CTCG_R2_subset.fq -R '@RG\tID:FLOWCELL1.LANE6\tSM:Z23337' | samtools view -Shu - | samtools sort - -o Z23337_sorted.bam
samtools index Z23337_sorted.bam
```

As you can see in the [`bwa` manual](http://bio-bwa.sourceforge.net/bwa.shtml), the `mem` flag tells `bwa` to align the reads to the reference genome (or, using the `bwa` jargon, find the coordinates of the reference genome that map each read) using the BWA_MEM algorithm. There are several additional options you could specify if you want.  For example, you could use the `-M` flag to set a limit on the number of mismatches between a read and the reference genome or the `-t` flag to set the number of threads for the mapping.  Because we are mapping data from one species to a reference genome from another, we will not do this.



`bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:XXX\tPL:illumina" reference.fa data.sai data.fq > data.sam | samtools view -bShu - > data_ref.bam`

The `samse` command of `bwa` tells the software that we are working with single-end data. The -r tells the software we wish to provide a readgroup (which is used later by the Genome Analysis Tookkit software). We then pipe the output using the vertical bar to the program `samtools`, which then prints the alignment (`view`) to an output file, The function of the `-bShu` flags can be found by typing `samtools view` at the command prompt.

For example, you could type this:

`bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:`SAMPLEID`\tPL:illumina" my_monkey_chromosome/chr`ZZZ`.fa samples/PF515_chr`ZZZ`.sai samples/PF515.fq | samtools view -bShu - > samples/PF515_chr`ZZZ`.bam`

Here you need to change the "ZZZ" to match the chromosome you are working with and the "SAMPLEID" to the sample you are working with (e.g. PF515). You also need to type the path and filenames of your reference genome, the .sai file, and the fastq file.

Now sort the `.bam` file:

`samtools sort samples/data.bam samples/data_sorted`

Here you do not need the `.bam` suffix after the `data_sorted` prefix; this suffix is added automatically by `samtools`. For example:

`samtools sort samples/PF515_chrX.bam samples/PF515_chrX_sorted`

Now please make an index for the bam file, which is a `.bai` file:

`samtools index samples/data_sorted.bam`

## Updates to bwa

The `mem` algorithm is an update to bwa that has a simpler pipeline for preparing the reference genome. For that, all you need to do is index the genome, and then type this:

`bwa mem -M -t 16 -r "@RG\tID:FLOWCELL1.LANE6\tSM:PF515\tPL:illumina" my_monkey_chromosome/chr9.fa samples/PF515.fq | samtools view -bSh - > samples/PF515.bam`

We have not used this algorithm for this example because the Stacks software we are using is not compatible with the output (although a recent update is).

## Practice Problem 4: Assessing coverage

Samtools can provide information on the number of reads for each position of the reference sequence for which there are data.  You can see this information by typing this:

`samtools depth XXX_sorted.bam`

Where `XXX` is the sample ID number.  If you want to know the average depth across all sites, you could type this:

`samtools depth XXX_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'`

Here, as previously, the vertical bar `|` is a "pipe" that sends the information from the command before it to the command after it.  So the data you generated from `samtools` will be parsed with the unix `awk` command.  This will add the values of the third column `$3` to a variable called `sum` and then at the end (`END`) print out the word `Average` followed by the quotient `sum/NR` where `NR` is a built in variable that keeps track of the number of records.  A good description of `awk` is [here](http://www.folkstalk.com/2011/12/good-examples-of-awk-command-in-unix.html).

## Practice Problem 5 (for home): De-multiplexing the complete dataset and mapping the data to your reference chromosome for one individual

Using the same pipeline we have just gone through, please do the following:
* demultiplex the complete dataset
* rename the resulting fastq files to match the sample names instead of the barcode sequences
* map the complete data from one individual (e.g. PF515) to your reference chromosome.

Now, using the idxstats options of [samtools](http://www.htslib.org/doc/samtools-0.1.19.html) please check how many reads mapped to your chromosome for your mapped data and how many failed to map.  Do you know why so many failed to map?

## OK, if this all went smoothly we are now ready to automate the alignments with a bash script.  Please click [here](https://github.com/evansbenj/BIO720/blob/master/3_Lecture_3_Automating_alignment_with_bash.md) to go to the next page.
