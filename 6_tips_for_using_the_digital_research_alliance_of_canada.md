# Tips for using ComputeCanada

Or you can go back to [Genotyping with bcftools](https://github.com/evansbenj/2024_BIO722/blob/master/5_genotyping_with_bcftools.md)

The [Digital Research Alliance of Canada](https://ccdb.alliancecan.ca/security/login) is an amazing resource!  To access this resource, your supervisor will need to register themself and you. Once you are registered, you have immediate access to hundreds of computer nodes, and hundreds of processing hours. If you need more, your supeprvisor can apply for more.

This resource uses a queueing system to run jobs.  The main commands I use are:
* `sbatch bashfile` where bashfile is a text file that has information that the queueing system needs to run a job
* `squeue -u username` where username is your login name
* `scancel jobid` where jobid is the identification number of a job you wish to terminate
* `module spider program_name` where program_name is the name of a program you want to use
* `module load dependencies program_name` where dependencies are the required tools that are needed to run your program (which is determined using the `module spider` command above

# Tips for setting up your bashfile

Here is an example of a typical bashfile that I would use to run a job
```sh
#!/bin/sh
#SBATCH --job-name=readgroups
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=readgroups.%J.out
#SBATCH --error=readgroups.%J.err
#SBATCH --account=def-XXX

# Here I use comments to remind me how to run the script. I am passing in a path of a directory
# sbatch ./2024_picard_add_read_groups.sh directory 

# Here I am loading the dependencies and software required for this job; this particular
# program has no dependencies but I do need to load a specific version which was
# determined by entering "module spider picard"
module load picard/2.23.3

# Now I am running the job in a loop. The ${1} refers to the directory I passed in with the sbatch command
for file in ${1}*_sorted.bam
do
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${file} O=${file}_rg.bam RGID=4 RGLB=$(basename $file) RGPL=ILLUMINA RGPU=$(basename $file) RGSM=$(basename $file)
done

```

The `--account=def-XXX` part needs to be updated with the username of your supervisor

There are several convenient features of this bash script. The job will run if you remember to load the module before you start the job, but it is easier (and more foolproof) to just load the module inside the script. Another nice feature is being able to run the script from any directory because you pass in the directory with your files in it. I usually keep all my sbatch scripts in one folder and then run the `sbatch` command from within the directory that has my files. The output files will be written to the directory from which you launch the `sbatch` command.

Some of the job parameters (such as `--time` and `--mem` require some trial and error (or knowledge that I do not have). 

The script also will create two files that have the jobID within the file name. For this script, these files will be called `readgroups.JOBID.err` and  `readgroups.JOBID.out`. The JOBID is a unique identifier that is assigned once you run the sbatch command. The nice thing about including the JOBID in your error and out files is that you can compare these files over multiple runs and assess whether you fixed a problem if there is one.
