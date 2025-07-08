def grab_samples(fastq_folder):
	'''
	Grabs sample names from an input fastq folder and returns them
	'''
	import os
	file_names=os.listdir(fastq_folder)
	samples=[]
	for file_name in file_names:
		samples.append(file_name.replace('.fastq.gz', ''))
	return samples

all_samples=grab_samples('fastq')

rule all:
	input:
		aligned_sam=expand('sam_files/{sample}_output.sam', sample=all_samples)

rule run_minimap:
	input:
		genome='Pf_3D7_genome/PlasmoDB-62_Pfalciparum3D7_Genome.fasta',
		sample='fastq/{sample}.fastq.gz'
	output:
		aligned_sam='sam_files/{sample}_output.sam'
	shell:
		'minimap2 -a -o {output.aligned_sam} {input.genome} {input.sample}'
