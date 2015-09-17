# Automating alignments to a reference genome with a bash script

(Or you can go back to [mapping reads to a reference genome](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/4_Mapping_reads_to_a_reference_genome.md)).

Now that you have seen how to align data from one individual to a reference genome, we can automate the alignment of all individuals to the reference genome using a bash script.  This is much better than going through all that stuff independently for each individual.  We can accomplish this by defining an `array` that contains the names of all of the individuals in the analysis, and then looping through this array and executing each of the commands for each individual.

To do this, we can use a bash script.  Let's first make sure we are still in the correct directory:

`pwd`

If you don't see `~/monkey`, please type this:

`cd ~/monkey`

Now we can make a bash script in this directory.  Here is an example that should run all of our analyses for each individual:

```
#!/bin/bash                                                                                                                  

path_to_bwa="/apps/bwa/0.7.12"
path_to_samtools="/apps/samtools/0.1.19"
path_to_data="."
path_to_chromosome="/home/datasets/2015_Ben_Evans/rhesus_chromosomes"
chromosome="chrXXX.fa"

individuals="PF515
PM561
PM565
PM566
PM567
PM582
PM584
PM592
PM602"

for each_individual in $individuals
do

echo ${each_individual}
    $path_to_bwa/bwa aln $path_to_chromosome/$chromosome $path_to_data/${each_individual}.fastq > $path_to_data/${each_individual}.sai
    $path_to_bwa/bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:${each_individual}.fastq\tPL:illumina" $path_to_chromosome/$chromosome $path_to_data/${each_individual}.sai $path_to_data/${each_individual}.fastq > $path_to_data/${each_individual}.sam
    $path_to_samtools/samtools view -bt $path_to_chromosome/$chromosome -o $path_to_data/${each_individual}.bam $path_to_data/${each_individual}.sam
    $path_to_samtools/samtools sort $path_to_data/${each_individual}.bam $path_to_data/${each_individual}_sorted
    $path_to_samtools/samtools index $path_to_data/${each_individual}_sorted.bam
done

```

In the beginning of the script 5 variables are defined that specify, respectively, the path for the bwa and samtools programs, the path to the data, the path to the reference chromosome, and the name of the chromosome you are working on.  Please copy this section and then make a new file called `alignment_commando` using emacs by typing this:

`emacs alignment_commando`

This should open up an emacs window.  You can then paste the bash script into this window.  Now use the arrow keys to scroll up to the line that says `chromosome="chrXXX.fa"` and change the part that says `chrXXX.fa` to correspond with whatever chromosome you are working on.  For example, if youa re working on chromosome 10, please change this to instead read `chr10.fa`.

Now type `Ctrl-x` and `Ctrl-s` to save the file.

Now we need to make the file executable, so type this:

`chmod +x alignment_commando`

And now we should be able to execute the file.  Type this:

`./alignment_commando`

## Problem 4

Imagine you were using this script to work with a complete genome alignment.  This could take some time and you want to have some idea of how much progress the script has made.  Can you please use the Unix `echo` command to the bash script to keep you informed about which command is executing.  You could (for example) ask the script to tell you when command 1, 2...5 is done and for which individual it has been completed.

Once this is done, we can try working with the program [Stacks](http://creskolab.uoregon.edu/stacks/manual/).  Stacks is actually a suite of programs will help us compile these data, calculate summary statistics, and also output the data to other software for further analysis.  

## Please click [here](https://github.com/evansbenj/Reduced-Representation-Workshop/blob/master/6_Using_Stacks_with_a_reference_genome.md) to go to the next page on using Stacks.
