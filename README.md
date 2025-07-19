# clair3-nanopore

This repository is being written collaboratively as part of a tutorial.
Eventually it will hold complete instructions for the following steps:

- concatenating reads with a custom function
- aligning reads with minimap2
- variant calling with clair3
- manipulating clair3 outputs to contain variant calls for variants of interest
- reformatting output vcf files for compatibility with SNPEff

The current version takes a (hardcoded) fastq folder and (hardcoded) genome
folder and runs minimap2 on all samples in parallel.

## temporary manual clair3

This is a temporary folder that contains example commands that a user can run in
order to run clair3 manually. This folder assumes that you have already run
minimap2 using the align_minimap.smk script.

The preparing_bam_reads.sh file includes example commands for indexing a genome
(the faidx command), converting a sam file to a bam format (the samtools view
command), sorting the sam file (the samtools sort command) and indexing the bam
file (the samtools index command).

The running_clair3.sh script includes example commands for running clair3 on an
example sample. This temporary folder will go away once we have incorporated
these commands into our main snakemake script.