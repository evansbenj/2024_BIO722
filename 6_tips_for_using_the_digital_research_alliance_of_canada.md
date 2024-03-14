Or you can go back to [Genotyping with bcftools](https://github.com/evansbenj/2024_BIO722/blob/master/5_genotyping_with_samtools_and_bcftools.md)

# Tips for using ComputeCanada

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

```
