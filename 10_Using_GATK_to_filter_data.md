# Using GATK to filter data

It is often the case, despite our efforts to generate high quality genotype calls, that some genotypes just don't make sense.  For example, we might observe a heterozygous genotype on the male specific portion of the Y chromosome or we might see some genotypes from a female on the Y chromosome. We can easily identify and screen out these sites (i.e. filter them) using `GATK`.
