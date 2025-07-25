configfile: 'align_minimap.yaml'
output_folder=config['output_folder']
fastq_folder=config['fastq_folder']

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

all_samples=grab_samples(fastq_folder)

rule all:
	input:
		aligned_sam=expand(output_folder+'/sam_files/{sample}_output.sam', sample=all_samples)

rule run_minimap:
	input:
		genome=config['genome_fasta'],
		sample=fastq_folder+'/{sample}.fastq.gz'
	output:
		aligned_sam=output_folder+'/sam_files/{sample}_output.sam'
	shell:
		'minimap2 -ax map-ont -o {output.aligned_sam} {input.genome} {input.sample}'
