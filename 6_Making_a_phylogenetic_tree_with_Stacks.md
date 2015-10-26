## Using `Stacks` to export a file for phylogenetic analysis.

(Or you can go back to using Stacks and POPBAM [here](https://github.com/evansbenj/BIO720/blob/master/5_more_on_Stacks.md)).

We can also use `Stacks` to export an input file for phylogenetic analysis.  Let's try this now.

From within your `Stacks_Results` directoyr, please use your favorite text editor to make a population map file called `population_map_all_different` in which all individuals were assigned to a different population.  This tab-delimited file looks like this:

```
PF515_sorted    PF515
PM561_sorted    PM561
PM565_sorted    PM565
PM566_sorted    PM566
PM567_sorted    PM567
PM582_sorted    PM582
PM584_sorted    PM584
PM592_sorted    PM592
PM602_sorted    PM602
```

Please make sure this file is tab delimited.  Then please type this command:

`/usr/local/stacks/bin/populations -P ./ -b 1 -r 1 -t 36 -M population_map_all_different --phylip --phylip_var`

Most of these flags were discussed previously. The `--phylip` flag combined with the `--phylip_var` flag tells `Stacks` to output sites that are variable between and within populations in `Phylip` format, which is the format of the input file for Joe Felsenstein's `Phylip` package. This can be easily modified for other programs, such as the `nexus` format. The `-M` flag tells `Stacks` to use the new population map file in which each individual is assigned to a different population.

This should generate a file called `batch_1.phylip`.  

Let's have a look at this file now.  Please type this:

`more batch_1.phylip`

You should be able to see a Phylip formatted file. You can press the space bar to scroll down the file and press `q` when you want to quit.  The first line has the number of taxa (9 in this case) followed by the number of characters of data.  The next line has a taxon name (`1`, which corresponds to the first sample PF515) followed by the sequence data.  Some of these data are regular nucleotides (A, C, G, or T) and others are [IUPAC symbols](http://www.bioinformatics.org/sms/iupac.html) that indicate heterozygous SNPs (Y for C/T, R for A/G, etc).  The lines after this give data for the next samples.

## Making a Quick phylogeny using `Phylip`

Now we can use a program called `Phylip` to make a quick phylogenetic tree.  You can find out about the programs available in Phylip [here](http://evolution.genetics.washington.edu/phylip/phylip.html). We will make a neighborjoining tree because it is very quick to do. The Phylip programs that do this for us are called `dnadist` and `neighbor`.  Please open up a window on your browser to look at the commands for these programs.  Note that for your research I recommend to instead use maximum likelihood or Bayesian methods to make phylogenetic trees, for example using software such as [MrBayes](http://mrbayes.sourceforge.net/), [BEAST](http://beast.bio.ed.ac.uk/), and [RaxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html).  We are using `Phylip` because it is quick.

To make a phylogenetic tree with `Phylip`, we first need to make a matrix of pairwise distances from our sequence data.  You can make a distance matrix by typing this:

`dnadist`

This will open up a menu of options.  You can enter `batch_1.phylip` for your input file.  I suggest naming your output file, which is the matrix of pairwise distances, something useful such as `batch_1.phylip.dist`.  You will use this as your input for the `Phylip` program called `neighbor`.

To run `neighbor` please type this (using the US spelling of neighbour):

`neighbor`

You will be prompted for a distance file and you will create two output files.  One is a not very useful graphical representation of the neighborjoining tree.  The second is a parenthetical notation of a phylogeny, including branch lengths.  Name this second file something useful such as `batch_1.phylip.tre`.

Now please download this tree file to your local computer using the scp command. I usually do this by opening up a new termina session and typing something like this:

`scp username@info.mcmaster.ca:path_to_the_file_I_want/filename directory_on_my_local_computer_where_I_want_the_file_to_go`

You can now view the tree you made using the [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) program. Please install this program if you haven't already.

# OK now let's use `Stacks` to make an input file for `Structure` analysis.  Please click [here](https://github.com/evansbenj/BIO720/blob/master/7_Stacks_and_Structure.md).

