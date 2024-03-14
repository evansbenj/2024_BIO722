# Quality control and trimming

(Or you can go back to the Readme page [here](https://github.com/evansbenj/BIO720/blob/master/0_README.md)).

## A quick note about genetic samples

Your statistical power, precision, and accuracy will depend on the quantity and quality of your data.  These factors, of course, depend on the starting material you use for sequencing.  Excitingly, a bunch of recent studies suggest that whole genome sequencing and reduced representation genome sequencing can be used on sub-optimal samples, for example from museum specimens, fecal extractions, and the like.  Nonetheless, it is in your interest to use as high quality DNA as possible.  
- For animal tissue, I recommend using ethanol (preferable for DNA) or RNAlater to preserve your tissues and I suggest chilling your samples (at -20 or -80 degrees) as soon as possible after collection.
- For DNA extraction, I've had success using Qiagen DNEasy extraction kits.  These kits can be used in any lab that has a heat block and a centrifuge.  When a small amount of starting material is being used I recommend reducing the volume of elution buffer you use in order to concentrate the DNA.  You can always dilute it later, and it is harder and less efficient to make a sample more concentrated.
- If you want to outsource the library preparation (a good idea if you can afford it), I recommend using an agarose gel to first roughly normalize the concentration of your samples and ensure that the quality is as good as possible (high molecular weight).  Different library preparation methods have different requirements in terms of the volume and concentration of gDNA.
- If you have many samples, I suggest extracting and comparing more than you plan to run so you can choose the best ones.

## How much does this cost and where can I do it?

I've recently had ddRADseq libraries made at [Ibis]([http://www.floragenex.com/](https://www.ibis.ulaval.ca/en/services-2/genomic-analysis-platform/)), which is located at Laval University in Quebec City. I've had sequencing done at [Genome Quebec](https://genomequebec.com/en/) in Montreal.  In 2024, the prices for a 96 RADseq sample run, including library preparation and normalization was ~CAD$4500. Sequencing with Illumina 150bp paired end reads costs around CAD$2400. I anticipate that the cost of these services will decline over the next few years and that more sequencing centers will offer this service.

## A quick note about Markdown and Github

This website is written in a markup language called [Markdown](https://en.wikipedia.org/wiki/Markdown) and hosted by [Github](www.github.com).  I've found both of these tools to be easy to learn and very useful.

# Introduction to Quality Control, De-multiplexing, and Trimming of Illumina Data

## Fasta and Fastq format

Illumina sequence data is provided in a text file that is in a format called `fastq`.  This is a modification of another format called `fasta` in which each sequence has a header that begins with a `>` sign.  This is followed by the sequences.  Here is an example:

```
>example_sequence_in_fasta_format
ATGCGCGCGCTAGGCTCGCGATCGGGGAGCGCGAGCTGAGCTAGCGCGATGCGCCCCGAC
```

The format of `fastq` files is similar to `fasta` except that quality scores are included.  Each sequence has four lines (instead of two for `fasta` files).  The first begins with `@` followed by information about the sequence.  The second line is the nucleotide sequence. The third line is a `+` which may be followed by the same information that followed the `@` sign in the first line.  The fourth line is the quality values.  For the Illumina data we will be working with, these values range from 0–41 and are represented by single characters.  More details about fastq format is available [here](http://en.wikipedia.org/wiki/FASTQ_format).  Here is an example of a sequence in fastq format:

```
@HWI-ST724:202:D127MACXX:6:2114:13665:74490
TGCAGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAAATTTTTGGGCACAAAGAACCACAGAAAAAAAATGAAAA
+HWI-ST724:202:D127MACXX:6:2114:13665:74490
AFHJJJFIJJJJIJJJJJHFDDDDDB0530&0)00&)0&05007BDD############################################
```

In this sequence the number signs indicate low quality reads at the end (right side) of the sequence.  

## Set up a directory on scratch and make symbolic links

Please login to info, connect to info2020 (ssh info2020) and then navigate to the scratch directory as follows:

`cd /scratch/USERNAME/`

Where USERNAME is your username (e.g. gradstud12). Now make a directory for you to work in:

`mkdir Froggies`

Next, please switch to that directory (`cd ZZZ`) and make another directory for the fastq files:
`mkdir fq'

Enter this directory (`cd fq`).

Now let's makes symbolic links the full dataset: 
```
ln -s /home/ben/2024_BIO722/2022_pygmaeus/fq/*fq.gz .
```

Let's see how many lines one of the files is:
```
zcat Z23337_CTCG_R1.fq.gz | wc -l
```
There is a lot of data! (The number of lines divided by 4 is the number of reads)

Make a subset of a forward and reverse read that you can use to map against a reference genome:
```
zcat Z23337_CTCG_R1.fq.gz | head -n 100000 > Z23337_CTCG_R1_subset.fq
zcat Z23337_CTCG_R2.fq.gz | head -n 100000 > Z23337_CTCG_R2_subset.fq
```
(Note here that the number of lines (`-n 100000`) must be divisible by 4 to prevent truncating the information for the last read.

Check out the last four lines of these two files:
```
tail -n4 Z23337_CTCG_R1_subset.fq
tail -n4 Z23337_CTCG_R2_subset.fq
```
OK, now we have the data set up for us to work with.

## Quality Control
Before we do anything with individual sequences, it is a good idea to survey the overall quality of the data.  We can do this with many free tools; for this class we will use a program called [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  To run this program please type this:

`/usr/local-centos6/fastqc/fastqc Z23337_CTCG_R1_subset.fq`

This should give you some feedback about the analysis as it runs and generate an `.html` file called `Z23337_CTCG_R1_subset_fastqc.html`.

Please download the `html` file to your local computer (later) and open it in a browser to view the quality metrics. You can do this by opening up a new terminal window (or Mobaxterm window) and typing something similar to this:
```
scp USERNAME@info.mcmaster.ca:PATH_TO_FILE/Z23337_CTCG_R1_subset_fastqc.html .
```
You can figure out the path by typing 'pwd'. (Note that you should leave out the `/2` in the beginning of your path.)

## Trimming

Let's use [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic) to remove poor quality positions in our data.

```
java -jar /usr/local-centos6/trimmomatic/trimmomatic-0.36.jar PE Z23337_CTCG_R1_subset.fq Z23337_CTCG_R2_subset.fq output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:/usr/local-centos6/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
```

## Practice problem 1 (for home): How many reads do we have for each individual?

As an exercise, please use the [`zgrep`](http://unixhelp.ed.ac.uk/CGI/man-cgi?zgrep) command to count how many reads we have for each of the .fq.gz files.  A hint is that using `zgrep`, you can count the number of times an identifier character for each sequence appears in each file for each individual.  Another hint is that you can get the manual for any `Unix` command by typing `man command`.  Which individual has the most reads?  Which has the least reads?  Can you think of a reason that some samples have lots of reads while others have less?  

How would your `grep` command differ for a `fasta` file compared to a `fastq` file?

## Practice problem 2: How did trimming affect the quality assay?

Please use FastQC to evaluate read quality of your trimmed sequences.  How do the trimmed reads differ from the untrimmed reads?

## OK, now we are ready to move on to mapping reads to a reference genome.  Please click [here](https://github.com/evansbenj/BIO720/blob/master/2_Lecture_2_reference_genomes_and_read_mapping.md).

