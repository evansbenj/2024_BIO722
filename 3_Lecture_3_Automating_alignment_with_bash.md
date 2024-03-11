# Automating alignments to a reference genome with a bash script

(Or you can go back to [mapping reads to a reference genome](https://github.com/evansbenj/BIO720/blob/master/2_Lecture_2_reference_genomes_and_read_mapping.md)).

Now that you have seen how to align data from one individual to a reference genome, we can automate the alignment of all individuals to the reference genome using a bash script. This is easier than going through all that stuff independently for each individual. We can accomplish this with a `bash` script that loops through the names of our fq files and executes each of the commands for each individual and names our bam files with information that is in the names of the fq files.

Below is an example `bash` script that can run all of our analyses for each individual.  Please use a text editor to make this program.  In the beginning of the script there is a "shebang" that tells the computer that this is a bash script. Then there is a loop (beginning with `for`) that asks the computer to read in all files from a directory that we passed in as an argument ($1) that have the suffix `_R1.fq.gz`. It then will output commands to the screen using echo. I have intentionally only output the commands to prevent you from accidentally running these commands. If you really wanted to run this, then just remove both `echo` commands.

```
#!/bin/bash                                                                    

for file in $1/*_R1.fq.gz ; do
    echo bwa mem ../reference/XENLA_10.1_genome.fa.gz $1/${file::-9}_R1.fq.gz $1/${file::-9}_R2.fq.gz -R "@RG\tID:FLOWCELL1.LANE6\tSM:${file::-9}" | samtools view -Shu - | samtools sort - -o ${file::-9}_sorted.bam
    echo samtools index ${file::-9}_sorted.bam
done

```


Now we need to make the file executable, so type this:

`chmod +x alignment_commando`

And now we should be able to execute the file.  Type this:

`./alignment_commando`

Here we needed to preceed the name of our `bash` script by `./` to tell the computer where to find our script (i.e. in the current working directory).

You should see some errors, but this is because the second command is expecting the output of the first command, but the `echo` prevented this file from being made.



##  Practice Problem 6 (for home): Modify this bash script to batch process other operations

Can you modify the bash script above to work with the subsetted fastq file you made?  Or with a different reference genome that you pass in as an argument?

## Practice Problem 7 (for home): Using a bash script to get coverage statistics from all individuals

Can you make a bash command to run trimmomatic on all of the paired-end data?

## Please click [here](https://github.com/evansbenj/BIO720/blob/master/4_Using_angsd_to_analyze_your_bam_files.md) to go to the next page on using Stacks.
