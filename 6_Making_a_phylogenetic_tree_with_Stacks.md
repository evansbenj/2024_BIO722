## Using `Stacks` to export a file for phylogenetic analysis.

(Or you can go back to Stacks and Structure [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/8_Stacks_and_Structure.md)).

Much in the same way that we used `Stacks` to export an input file for the program `Structure`, we can also use `Stacks` to export an input file for phylogenetic analysis.  Let's try this now.

Because this (again) takes a little while, Ben did it in advance using the following commands (You do not need to type these).  First he went to the results directory from the complete dataset:

`cd /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results`

Then he made a population map file called `population_map_all_different` in which all individuals were assigned to a different population.  This tab-delimited file looks like this:

```
PF515_sorted    population_1
PM561_sorted    population_2
PM565_sorted    population_3
PM566_sorted    population_4
PM567_sorted    population_5
PM582_sorted    population_6
PM584_sorted    population_7
PM592_sorted    population_8
PM602_sorted    population_9
```

...and then he typed this command:

`/apps/stacks/1.29/bin/populations -P ./ -b 1 -r 1 -t 36 -M population_map_all_different --phylip --phylip_var`

Most of these flags were discussed previously.  The `--phylip` flag combined with the `--phylip_var` flag tells `Stacks` to output sites that are variable between and within populations in `Phylip` format, which is the format of the input file for Joe Felsenstein's `Phylip` package.  This can be easily modified for other programs, such as the `nexus` format.  The `-M` flag tells `Stacks` to use the new population map file in which each individual is assigned to a different population.

This generated a file called `batch_1.phylip.tsv`.  

Please copy this file to your home directory like this:

`scp /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results/batch_1.phylip ~/monkey/Stacks_Results`

Let's have a look at this file now.  Please type this:

`more batch_1.phylip`

You should be able to see a Phylip formatted file. You can press the space bar to scroll down the file and press `q` when you want to quit.  The first line has the number of taxa (9 in this case) followed by the number of characters of data (70484 in this case).  The next line has a taxon name (`1`, which corresponds to the first sample PF515) followed by the sequence data.  Some of these data are regular nucleotides (A, C, G, or T) and others are [IUPAC symbols](http://www.bioinformatics.org/sms/iupac.html) that indicate heterozygous SNPs (Y for C/T, R for A/G, etc).  The lines after this give data for the next samples.

## Making a Quick phylogeny using `Phylip`

Now we can use a program called `Phylip` to make a quick phylogenetic tree.  You can find out about the programs available in Phylip [here](http://evolution.genetics.washington.edu/phylip/phylip.html). We will make a neighborjoining tree because it is very quick to do. The Phylip programs that do this for us are called `dnadist` and `neighbor`.  Please open up a window on your browser to look at the commands for these programs.  Note that for your research I recommend to instead use maximum likelihood or Bayesian methods to make phylogenetic trees, for example using software such as [MrBayes](http://mrbayes.sourceforge.net/), [BEAST](http://beast.bio.ed.ac.uk/), and [RaxML](http://sco.h-its.org/exelixis/web/software/raxml/index.html).  We are using `Phylip` because it is quick.

To make a phylogenetic tree with `Phylip`, we first need to make a matrix of pairwise distances from our sequence data.  You can make a distance matrix by typing this:

`/apps/PHYLIP/3.696/exe/dnadist`

This will open up a menu of options.  You can enter `batch_1.phylip` for your input file.  I suggest naming your output file, which is the matrix of pairwise distances, something useful such as `batch_1.phylip.dist`.  You will use this as your input for the `Phylip` program called `neighbor`.

To run `neighbor` please type this:

`/apps/PHYLIP/3.696/exe/neighbor`

You will be prompted for a distance file and you will create two output files.  One is a not very useful graphical representation of the neighborjoining tree.  The second is a parenthetical notation of a phylogeny, including branch lengths.  Name this second file something useful such as `batch_1.phylip.tre`.

You can copy it to your local computer by opening up another `PuTTy` session and typing this:

`scp username@caf-hpc1.sun.ac.za:~/monkey/Stacks_Results/batch_1.phylip.tre .`

You can now view the tree you made using the [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) program that should be available on your local computer.

## And this is just the beginning

We covered some simple examples aimed at getting you started.  Not surprisingly there is much more to learn and each individual project will require some customization. 
- To further increase your analytical flexibility, a simple first step would be to learn a simple computing language that can allow you to parse files and do simple calculations and iterative tasks.  I use `Perl` for this and `Python` is also a great option.  
- Another obvious resource is the internet.  If you want to do something or have a question about a software, chances are that someone else has already had the same goal or question.  
- If you are really feeling crazy, you can even read program manuals and README files that come with the software.



