# BIO722
This is the BJE portion of the BIO722 graduate course in Advanced Bioinformatics in 2024.  Perhaps overly optimistically, this section will be squeezed into one 2-hour interactive lecture where the class will explore genomic data, including:
* trimming fastq files
* preparing a reference genome
* mapping data to a reference genome
* doing a genome wide association study (GWAS)
* genotyping
* some discussion about genotype filtering

 

## Background on Reduced Representation Genome Sequencing
Many interesting organisms have big genomes, making complete genome sequencing infeasible, especially for multiple individuals.  A relatively cost-efficient solution has been developed recently called "reduced representation genome sequencing" (RRGS).  This approach enables deep sequencing of multiple genetic samples from the same (homologous) genomic regions.  It takes advantage of next generation (Illumina) sequencing technology and requires relatively simple laboratory preparation (DNA extraction), which can be accomplished with a centrifuge and a heat block.  Other laboratory steps (library construction) can be outsourced or done in house, depending on the equipment and funds that are available. This approach has many applications, including phylogenomics, population genomics, linkage mapping, and analysis of gene flow.

## Goals
The goal of this section of the course is to introduce students to some basic aspects of bioinformatic analysis of genomic data. Ideally this  will provide a sufficient level of exposure to students that they will be able to figure out how to learn more and generate and analyze their own datasets. Many of the approaches to analyze RRGS data are identical to WGS data (e.g., demultiplexing, trimming, mapping to a reference, read filtering, genotyping, etc).

## Interactive session 
To illustrate how to work with and map genomic data to a reference genome, we will first work with a RRGS dataset from the Bouchia clawed frog (*Xenopus pygmaeus*).  This frog occurs in central Africa, mostly in the Democratic Republic of the Congo. The data come from a captive-reared family where brothers and sisters were raised through metamorphosis, euthanized, and dissected to determine their sex. The goal of the project is to determine where the sex determining region of this species is. Using reduced representation genome sequencing, we will map the sex determining region of this species to a chromosome-scale reference genome from a closely related species (*Xenopus laevis*), and then explore what it can tell us about how important things evolve in general.  
 
## OK, lets begin by clicking [here](https://github.com/evansbenj/BIO720/blob/master/1_Lecture_1.md).
