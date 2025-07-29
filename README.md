# clair3-nanopore

The goal of this pipeline is to convert raw Nanopore reads into tables that list
the number of reads that support different amino acid mutations in each sample.
These AA tables are suitable for calculating prevalences, as well as for
graphing and visualizing mutations.

The visualization part is handled by another repo, here:
(https://github.com/simkin-bioinformatics/seekdeep_amplicon_visualization)

This clair3-nanopore repository is being written collaboratively as part of a
tutorial. Eventually it will automatically perform the following steps:

- concatenating reads with a custom function
- aligning reads with minimap2
- variant calling with clair3
- manipulating clair3 outputs to contain variant calls for variants of interest
- reformatting output vcf files for compatibility with SNPEff

## Current progress:

- concatenated reads are aligned with minimap2
- alignments are converted to bam, sorted, and indexed
- genome gets indexed with samtools faidx
- aligned reads are variant called with clair3 to produce VCF outputs
- a custom script copies gff and genome files for snpEff
- snpEff builds a database from gff and genome files
- snpEff annotates VCF files

## Next Steps (not implemented yet):

- snpEff VCF gets converted to AA tables with custom script
- add ability to correctly parse "targeted" mutation depths
- incorporate snpEff converter into snakemake pipeline
- add ability to concatenate raw nanopore reads (near beginning of pipeline)
- more detailed yaml file instructions

## How to run the pipeline:

1. Clone this github repository.
2. Obtain and install a copy of mamba or micromamba (you can google this).
3. Create a new conda environment, e.g. with mamba env create -n clair3.
4. Activate your new environment, and install the following conda packages in
   your environment: clair3, snakemake, minimap2, bcftools, samtools.
5. Install snpEff (for me a quick google search revealed that this involves wget
   of snpEff and unzipping a folder).
6. cd into the github repository folder, edit the
   nanopore_variant_annotation.yaml file using the instructions provided in the
   yaml file, and run the program with this command:

`shell
snakemake -s nanopore_variant_annotation.smk --cores 4
`