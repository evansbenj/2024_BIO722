# Using a whitelists and blacklists to explore polymorphism inside and outside of genes.

## Downloading an annotation of rhesus genes

One advantage to working with an annotated genome sequence is that we can easily access information about gene locations. For the rhesus genome we are working with, we can download the coordinates of known genes [here](http://genome.ucsc.edu/cgi-bin/hgTables). Please visit this page now.

To download the annotations for rhemac2, please select "Mammal" from the clade menue, "Rhesus" from the genome menu, and "Jan. 2006 (MGSC Merged 1.0/rhemac2)" from the assembly menu.  In the group menu please select "Genes and Gene predictions" and for the track menu please select "RefSeq Genes".  In the table menu, please select "refGene" and in the region option please click the "genome" button.  Now change the output format to "BED - browser extensible data" and enter in the output file a name for your annotation file such as "rhemac2.bed". Once you have done all this, please click the "get output" button. Please upload this to your info account.

You can check out the contents of this file by typing:

`more rhemac2.bed`

This [page](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) provides information about many genomic formats, including `bed` files. Please check it out.  As you can see, the first three columns are the ones we care about: these tell us the chromosome number, start and end coordinates of genes.  Please note that the numbering format for bed files is somewhat weird; the start coordinate is based on a numbering system beginning with zero, not one. And the end coordinate refers to the base after the end of the feature. As stated on this page, the first 100 base pairs of a chromosome would have bed coordinates of 0 to 100, and would span bases that are numbered 0 to 99.

## Making a whitelist from the bed file for use with `Stacks`

I (hopefully) wrote a perl script that hopefully will generate a list of RAD loci from the `batch_1.catalog.tags.tsv` file that fall within the intervals provided in the `rhemac2.bed` file. We can then use this to generate polymorphism statistics inside and outside of genes with the `populations` module of `Stacks`.

