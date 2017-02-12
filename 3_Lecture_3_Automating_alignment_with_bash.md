# Automating alignments to a reference genome with a bash script

(Or you can go back to [mapping reads to a reference genome](https://github.com/evansbenj/BIO720/blob/master/2_Lecture_2_reference_genomes_and_read_mapping.md)).

Now that you have seen how to align data from one individual to a reference genome, we can automate the alignment of all individuals to the reference genome using a bash script. This is easier than going through all that stuff independently for each individual. We can accomplish this with a `bash` script by defining an `array` that contains the names of all of the individuals in the analysis, and then looping through this array and executing each of the commands for each individual.

Below is an example `bash` script that should run all of our analyses for each individual.  Please use a text editor to make this program.  In the beginning of the script 5 variables are defined that specify, respectively, the path for the bwa and samtools programs, the path to the data, the path to the reference chromosome, and the name of the chromosome you are working on. You will need to modify the variables somewhat to match the chromosome you are working on and the directory. For example you should use the arrow keys to scroll up to the line that says `chromosome="chrXXX.fa"` and change the part that says `chrXXX.fa` to correspond with whatever chromosome you are working on.  For example, if youa re working on chromosome 10, please change this to instead read `chr10.fa`. Also, in the `path_to_chromosome` variable, you will need to change the part that says `YYY` to match your home directory name.

```
#!/bin/bash                                                                                                                  

path_to_data="/home/YYY/my_monkey_data"
path_to_chromosome="/home/YYY/my_monkey_chromosome/"
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
    bwa aln $path_to_chromosome/$chromosome $path_to_data/${each_individual}_subset.fastq > $path_to_data/${each_individual}.sai
    bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:${each_individual}_subset.fastq\tPL:illumina" $path_to_chromosome/$chromosome $path_to_data/${each_individual}.sai $path_to_data/${each_individual}_subset.fastq > $path_to_data/${each_individual}.sam
    samtools view -bt $path_to_chromosome/$chromosome -o $path_to_data/${each_individual}.bam $path_to_data/${each_individual}.sam
    samtools sort $path_to_data/${each_individual}.bam $path_to_data/${each_individual}_sorted
    samtools index $path_to_data/${each_individual}_sorted.bam
done

```



Now we need to make the file executable, so type this:

`chmod +x alignment_commando`

And now we should be able to execute the file.  Type this:

`./alignment_commando`

Here we needed to preceed the name of our `bash` script by `./` to tell the computer where to find our script (i.e. in the current working directory).

Now please check whether this worked by checking out the file size of the files in your `~/my_monkey_data` directory like this:

`ls -l ~/my_monkey_data`

You should see lots of file including, for each sample a `_sorted.bam` file that is not of file size zero.

##  Practice Problem 6: Using a bash script to mapping the complete data from all individuals

Imagine you were using the script above to work with a complete genome alignment.  This could take some time and you want to have some idea of how much progress the script has made.  Can you please use the Unix `echo` command to the bash script to keep you informed about which command is executing.  You could (for example) ask the script to tell you when command 1, 2...5 is done and for which individual it has been completed.

Once this is done, we can try working with the program [Stacks](http://creskolab.uoregon.edu/stacks/manual/).  Stacks is actually a suite of programs will help us compile these data, calculate summary statistics, and also output the data to other software for further analysis.  

## Practice Problem 7: Using a bash script to get coverage statistics from all individuals

Please write a bash script that will provide the coverage statistic for all bam files in a folder that have the suffix "_sorted.bam".  You can (and should) use the internet for tips.

## Please click [here](https://github.com/evansbenj/BIO720/blob/master/4_Using_Stacks_to_analyze_your_bam_files.md) to go to the next page on using Stacks.
