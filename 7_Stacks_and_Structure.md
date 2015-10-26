# Using `Stacks` to export an input file for `Structure`

(Or you can go back to using `Stacks` to calculate summary statistics [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/7_More_on_Stacks.md)).

[`Structure`](http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/html/structure.html) is a software that attempts to assign individuals to *k* populations in such a way as to minimize Hardy-Weinberg and linkage disequilibrium.  We run `Structure` by specifying multiple values of *k* and then seeing which value(s) maximuze the likelihood of the data given the model of population structure. We can generate an input file for this program using `Stacks`. 

Because this takes a while, Ben previously generated this file for us to work with.  He did this using the following commands (please don't type these).  First he made sure he was in the Stacks_Results directory that is from the complete dataset:

`cd /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results`

The program `Structure` can not handle all of our data from the complete dataset, so Ben selected 1000 loci randomly to analyze like this:

`shuf -n 1000 /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results/batch_1.catalog.tags.tsv | awk '{print $3}'  > 1000_randoms`

This uses another `unix` command called `shuf`.  This tells the computer to print a randomly selected value from column 3 1000 times to a file called `1000_randoms`.

Then Ben was ready to generate an input file for `Structure`.  Ben typed this command:

`/apps/stacks/1.29/bin/populations -P ./ -b 1 -r 1 -t 36 --structure --write_single_snp -W ./1000_randoms -B chrX_list`

This command directed the `populations` module of `Stacks` to output a single snp (the `--write_single_snp` flag) from tags specified by the `1000_randoms` file (the `-W` tag) but not to include any SNPs from chromosome X (the `-B` flag).

This generated a file called `batch_1.structure.tsv` which can be used as an input file for the program `Structure`.  Please copy this to your home directory like this:

`scp /home/datasets/2015_Ben_Evans/complete_data/monkey/Stacks_Results/batch_1.structure.tsv ~/monkey/Stacks_Results`

## Running Structure

Before we can run `Structure` with out data, we have some clerical issues to take care of.  First, please make sure you are in your home directory by typing this:

`cd ~/monkey/Stacks_Results`

Now lets delete the first row from our input file  - this row has a comment that we don't need.  Please type this:

`tail -n+2 batch_1.structure.tsv > simple_structure.tsv`

This uses the Unix command `tail` to feed all lines beginning with the second line to a new file called `simple_structure.tsv`.

Because we required that there be no missing data from the SNPs in our analysis, we can get a differing number of SNPs in our input file depending on which sites were randomly chosen. So we need to count how many columns we have in this file.  Please type this:

`head -n1 simple_structure.tsv |  sed 's/\t/\n/g' | wc -l`

This should output the number of tab spaces in the first column of the file `simple_structure.tsv`.  Because the first column begins with a tab space, the actual number of loci is equal to this number **minus one**.

Now we are ready to run `Structure`.  Please type this command:

`/apps/Structure/bin/structure -m /apps/Structure/sample/mainparams -e /apps/Structure/sample/extraparams -K 3 -L 268 -N 9 -i simple_structure.tsv -o output_K_3`

This tells the system to execute the `Structure` program and it specifies the paths to two files (`mainparams` and `extraparams`) that are used in the analysis.  It then has flags for the number of populations (`-K`), the number of loci (`-L`; based on the number we got above from the `head` command), the number of individuals (`-N`; our study has 9 individuals), the input file (`-i`) and an output file (`-o`).

 Assuming the command executes without error, you can check out the results in the file `output_K_3` like this:
 
 `more output_K_3_f`
 
 You could scroll down to the population assignments, which should look something like this:
 ```
 Inferred ancestry of individuals:
        Label (%Miss) Pop:  Inferred clusters
  1 PF515_sorte    (0)    1 :  0.687 0.150 0.164 
  2 PM561_sorte    (0)    1 :  0.583 0.194 0.222 
  3 PM565_sorte    (0)    1 :  0.635 0.198 0.167 
  4 PM566_sorte    (0)    1 :  0.665 0.171 0.164 
  5 PM567_sorte    (0)    1 :  0.591 0.205 0.204 
  6 PM582_sorte    (0)    1 :  0.674 0.149 0.177 
  7 PM584_sorte    (0)    1 :  0.613 0.184 0.203 
  8 PM592_sorte    (0)    1 :  0.591 0.200 0.209 
  9 PM602_sorte    (0)    1 :  0.521 0.226 0.253 
```

This tells us, for each individual, what the probability that that individual is assigned to each one of *k* populations (three in this case).

## A quick example of plotting with `R`

To use `R` on this system, we need to load a module.  Please type this:

`module load app/R`

Now we can plot this by making a file and pasting these data in this file:

`emacs assignments`

If you now paste the data beginning with the line `1 PF515_sorte...` and save it (`Ctrl-X` and then `Ctrl-S`) and exit (`Ctrl-X` and then `Ctrl-C`) we can easily plot the data with R.  To do this type: `R`.  This should open up the `R` environment.  Now import the data:

`>` `dat<-read.table(file="assignments")`

and make a pdf..

`>`  `pdf("temp.pdf", height=3, width=6)`

and make a barplot

`>` `barplot(as.matrix(t(dat[,6:8])), names.arg=dat$V2, col=rainbow(3), las=2, cex.names=0.65, cex.axis=0.65)`

This command tells R to make a boxplot using columns 6 thru 8 of the table called `dat`.  We have transformed these data for this plot using `t()`.

and now exit the `R` environment:

`>` `q()`

We can now dowload and view this file.  To download it, please open a new session and type this:

`scp username@caf-hpc1.sun.ac.za:~/monkey/Stacks_Results/temp.pdf .`

where `username` is your username.  You should be prompted for a password and then your file will download to your local computer.  You should be able to find this file, and then double click on the file and view it.

## Problem 6

OK, now please generate a plot from Structure with *k* equal to 2 and then another plot with *k* equal to 4.  When can you conclude by comparing these plots?

## OK, now let's make a [quick phylogeny using the RADseq data](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/9_Stacks_and_Phylogenies.md).
