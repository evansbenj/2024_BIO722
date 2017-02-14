# Calculating pairwise nucleotide diversity with Stacks

(Or you can go back to [Using Stacks with a Reference Genome](https://github.com/evansbenj/BIO720/blob/master/4_Using_Stacks_to_analyze_your_bam_files.md)).

If we want to calculate pairwise nucleotide diversity on the chromosome to which we mapped our data, we can use the `populations` module of `Stacks` as follows (from within the `Stacks_Results` directory):

`/usr/local/stacks/bin/populations -P ./ -b 1 -r 1.0 -t 36`

As detailed in the [manual](http://catchenlab.life.illinois.edu/stacks/comp/populations.php), this command directs the `populations` module of `Stacks` to write results to a directory specified by the `-P` flag. The `-r` flag says we want to only print data where 100% of the individuals have a genotype.  The `-b` and `-t` flags specify, respectively, that `populations` should focus on batch_ID number 1 (you can work with special IDs if you need to but this is beyond the scope of this workshop) and that `populations` should use 36 threads to do the calculations.  

From within the `Stacks_Results` directory, we can view the pairwise nucleotide diversity statistic ($$\pi$$) in a file called `batch_1.sumstats_summary.tsv` like this:

`more batch_1.sumstats_summary.tsv`

As indicated in the header, $$\pi$$ should be listed in the sixth to last column on the bottom row. For autosomes this statistic is accurate but there are problems with this calculation using only the X chromosome because some of the individuals are males; Ben will discuss these problems and a workaround he did with a perl script. The estimate of $$\pi$$ he recovered from the Stacks analysis of the bam files is 0.00108 substutitons per site.

Now lets take some time for you to all analyze your bam files. Please let Ben know what your calculation of pi is for each of your chromosomes.

# Calculating divergence to the reference genome using `POPBAM`

As some of you may realize, the mutation rates in males and females are not equivalent. In many species, including all primates, more cell divisions occur to make sperm than eggs, and for this reason (we think) the mutation rate is higher in males than in females.  Because X chromosomes spend 2/3rds of their time in females, we expect lower polymorphism on the X by virtue of this.  We can control for differences in mutation rate very simply – by taking the ratio of pi to divergence in the X and the aDNA. We calculated pi above, now we need divergence.

[`POPBAM`](http://popbam.sourceforge.net/) is a program that is useful for calculating several population genetic statistics. Some of the functionality overlaps with that of `Stacks`.  We will use `POPBAM` to calculate the level of divergence of one of our samples to the reference genome.

First we need to add some additional information to the readgroup header of our `bam` file.  From within the `samples` folder, please type this:

```
samtools view -H ZZZ.bam > header.sam
perl -pi.old -e 's{PL:illumina}{PL:illumina\tPO:POP1}g' header.sam
samtools reheader header.sam XXX.bam > XXX_new.bam
samtools index XXX_new.bam
```

The first line uses samtools to make a text file called `header.sam` that contains the header information for the file called `XXX.bam`.  You should use one of your sorted bam files for this. The second line uses `Perl` to search and replace text in the header.sam file.  Basically this adds text to the readgroup portion of this file. Thie third line uses samtools to change the header of our file. And the fourth line makes a new index file for our new bam file.

Now we can get divergence information using `POPBAM` (from within the `samples` directory) like this:

`popbam diverge -o 0 XXX_new.bam chrX -f ../my_monkey_chromosome/chrXXX.fa chrXXX > divergence.txt`

Don't forget to specify the chromosome that you are surveying at the end (`chrXXX`). This will write the output of `POPBAM` to a text file called `divergence.txt`.  The divergence Ben got for PF515 was 0.00524.

Ben will now collect the data you collected, summarize it, compare it to the diversity and divergence estimates from the X, and discuss what implications it has for social system evolution in our monkeys.

# Using a whitelists and blacklists to explore polymorphism inside and outside of genes

One concern about these statistics is that genomic regions that contain genes may be subject to effects of natural selection. To test this, we can use annotation information for our reference genome to identify and exclude RADtags that are located within genes. To do this we will download annotations for the rhemac2 genome and use a perl script to generate a "blacklist" of genomic regions to exclude from our analysis.

## Downloading an annotation of rhesus genes

One advantage to working with an annotated genome sequence is that we can easily access information about gene locations. For the rhesus genome we are working with, we can download the coordinates of known genes [here](http://genome.ucsc.edu/cgi-bin/hgTables). Please visit this page now.

If you wanted to download the annotations for rhemac2, you could select "Mammal" from the clade menue, "Rhesus" from the genome menu, and "Jan. 2006 (MGSC Merged 1.0/rhemac2)" from the assembly menu.  Then, in the group menu you would select "Genes and Gene predictions" and for the track menu select "RefSeq Genes".  In the table menu, you would select "refGene" and in the region option click the "genome" button.  Now change the output format to "BED - browser extensible data" and enter in the output file a name for your annotation file such as "rhemac2.bed". Once you have done all this, you would click the "get output" button. If you want to download this and upload it to your account, go ahead.  Alternatively, within the `my_monkey_chromosome` directory, just make a symbolic link to my version of this file like this:

`ln -s /home/ben/2015_BIO720/rhemac2.bed`

You can check out the contents of this file by typing:

`more rhemac2.bed`

This [page](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) provides information about many genomic formats, including `bed` files. Please check it out.  As you can see, the first three columns are the ones we care about: these tell us the chromosome number, start and end coordinates of genes.  Please note that the numbering format for bed files is somewhat weird; the start coordinate is based on a numbering system beginning with zero, not one. And the end coordinate refers to the base after the end of the feature. As stated on this page, the first 100 base pairs of a chromosome would have bed coordinates of 0 to 100, and would span bases that are numbered 0 to 99.

## Let's work with the full dataset mapped to one chromosome

To save time, I mapped the full data to chromosome 9 and sorted and indexed the bam files. Please make a symbolic link to these files like this:
```
ln -s /1/scratch/ben/PF515_chr9_sorted.bam
ln -s /1/scratch/ben/PM561_chr9_sorted.bam
ln -s /1/scratch/ben/PM565_chr9_sorted.bam
ln -s /1/scratch/ben/PM566_chr9_sorted.bam
ln -s /1/scratch/ben/PM567_chr9_sorted.bam
ln -s /1/scratch/ben/PM582_chr9_sorted.bam
ln -s /1/scratch/ben/PM584_chr9_sorted.bam
ln -s /1/scratch/ben/PM592_chr9_sorted.bam
ln -s /1/scratch/ben/PM602_chr9_sorted.bam
```
Then, please make a new Stacks directory and execute the Stacks pipeline for these data as described above.


## Making a whitelist from the bed file for use with `Stacks`

The `populations` module of `Stacks` can be fed a list of RAD loci to include (a whitelist) or exclude (a blacklist) from the analysis. The RADloci are listed in the file called `batch_1.catalog.tags.tsv`. Please unzip this file like this:

`gunzip batch_1.catalog.tags.tsv.gz`

and now check out the file like this:

`more batch_1.catalog.tags.tsv`

Each row corresponds to one RADlocus and the numbers in the third column are the ones we need to make a whitelist or a blacklist. The locations of the RADloci are indicated in columns 4 (the chromosome) and column 5 (the coordinate). Our task for each RADlocus is to check the annotation file (rhemac2.bed) to see if the RADlocus is inside a gene.

Fortunately we are not going to do this by hand. I wrote a perl script that will generate a list of RAD loci from the `batch_1.catalog.tags.tsv` file that fall within the intervals provided in the `rhemac2.bed` file. We can then use this to generate polymorphism statistics inside and outside of genes with the `populations` module of `Stacks`.

Here it is; please use your favorite text editor to make this file. Please change the permissions to allow it to be executable (`chmod +x filename`) and execute it as directed in the comment section. Ben will briefly explain what each part of the script does.

``` perl
#!/usr/bin/perl 
use strict;

# This program generates a whitelist for the populations module of Stacks
# given a catalog file (e.g. batch_1.catalog.tags.tsv) and a bed formated annotation
# file (e.g. rhemac2.bed)

# To execute this program please type this:
# make_my_whitelist.pl arg1 arg2 arg3
# where arg1 and arg2 are catalog and bed file names respectively and
# arg3 is an output filename
# for example:
# make_my_whitelist.pl batch_1.catalog.tags.tsv rhemac2.bed gene.whitelist

# first declare filename variables
my $inputfile = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $outputfile = $ARGV[2];

#### Prepare the input files 
unless (open DATAINPUT, $inputfile) {
	print "Can not find the catalog input file!\n";
	exit;
}

unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the annotation bed file!\n";
	exit;
}

# check the output file
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile, please provide a name!\n";
	exit;
}
print "Creating output file: $outputfile\n";

# define some variables
my %RADhash; # This is a data structure called a hash that will have format $RADhash{chromosome}[position] with value "locusnumber"
my @line; # this is a datastructure called an array
my $n; # this is a data structure called a string
my $rad_positions; # another string

# Now we are going to loop through each line of the catelog file
while ( my $line = <DATAINPUT>) {
	# split the line by tabs and load up an array called @line
	@line = split("\t",$line);
	# load the location and locus information into the RADhash
	$RADhash{$line[3]}{$line[4]}=$line[2];
} 
# Now we are going to loop through each line of the annotation file
while ( my $line = <DATAINPUT2>) {
	# split the line by tabs and load up an array called @line
	@line = split("\t",$line);
	# for each line, check if there is a $RADhash entry that lies between $line[1] and $line[2]
	for $rad_positions ( keys %{ $RADhash{$line[0]} } ){
		print $rad_positions,"\n";
		if (($rad_positions > $line[1])&&($rad_positions < $line[2])){
			# if yes, print the value of ($RADhash{$line[0]}{$rad_positions} to the outfile
			print OUTFILE $RADhash{$line[0]}{$rad_positions},"\n";
		}
	}
}

```
From within the `Stacks_Results` directory, I executed this script like this: 
`./make_my_whitelist.pl batch_1.catalog.tags.tsv ../../my_monkey_chromosome/rhemac2.bed gene.whitelist`

Assuming that went smoothly, you can now compare pairwise nucleotide diversity inside and outside of genes using these commands:

From within the `Stacks_Results` directory, you can calculate nucleotide diversity inside of genes using this `whitelist`, like this:

`/usr/local/stacks/bin/populations -P ./ -b 1 -r 1.0 -t 36 -W gene.whitelist`

Now check the pairwise diverisity within genes like this:

`more batch_1.sumstats_summary.tsv`

For nucleotide diversity outside of genes, use a `blacklist`:

`/usr/local/stacks/bin/populations -P ./ -b 1 -r 1.0 -t 36 -B gene.whitelist`

Now check the pairwise diverisity within genes like this:

`more batch_1.sumstats_summary.tsv`

Which was higher?
# Now let's use `Stacks` to make a phylogenetic tree [here](https://github.com/evansbenj/BIO720/blob/master/6_Making_a_phylogenetic_tree_with_Stacks.md).
