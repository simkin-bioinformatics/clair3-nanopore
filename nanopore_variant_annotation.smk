configfile: 'nanopore_variant_annotation.yaml'
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
sample_subset=['PM2', 'JS-2', 'JS-0-5-7']

rule all:
	input:
		#aligned_sam=expand(output_folder+'/sam_files/{sample}_output.sam', sample=all_samples),
		#indexed_bam=expand(output_folder+'/sorted_bam_files/{sample}_output_sorted.bam.bai', sample=all_samples),
		clair3_gvcf=expand(output_folder+'/clair3_samples/{sample}/merge_output.gvcf.gz', sample=sample_subset),
		copied_snakefile=output_folder+'/snakemake_params/nanopore_variant_annotation.smk'

rule copy_files:
	input:
		original_snakefile='nanopore_variant_annotation.smk',
		original_yaml='nanopore_variant_annotation.yaml'
	output:
		copied_snakefile=output_folder+'/snakemake_params/nanopore_variant_annotation.smk',
		copied_yaml=output_folder+'/snakemake_params/nanopore_variant_annotation.yaml'
	shell:
		'''
		cp {input.original_snakefile} {output.copied_snakefile}
		cp {input.original_yaml} {output.copied_yaml}
		'''

rule run_minimap:
	input:
		genome=config['genome_fasta'],
		sample=fastq_folder+'/{sample}.fastq.gz'
	output:
		aligned_sam=output_folder+'/sam_files/{sample}_output.sam'
	shell:
		'minimap2 -ax map-ont -o {output.aligned_sam} {input.genome} {input.sample}'

rule index_genome:
	input:
		genome=config['genome_fasta']
	output:
		genome_index=config['genome_fasta']+'.fai'
	shell:
		'samtools faidx {input.genome}'

rule sam_to_bam:
	input:
		aligned_sam=output_folder+'/sam_files/{sample}_output.sam'
	output:
		unsorted_bam=output_folder+'/unsorted_bam_files/{sample}_output.bam'
	shell:
		'samtools view -b -o {output.unsorted_bam} {input.aligned_sam}'

rule sort_bam:
	input:
		unsorted_bam=output_folder+'/unsorted_bam_files/{sample}_output.bam'
	output:
		sorted_bam=output_folder+'/sorted_bam_files/{sample}_output_sorted.bam'
	shell:
		'samtools sort -o {output.sorted_bam} {input.unsorted_bam}'

rule index_bam:
	input:
		sorted_bam=output_folder+'/sorted_bam_files/{sample}_output_sorted.bam'
	output:
		indexed_bam=output_folder+'/sorted_bam_files/{sample}_output_sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'

rule run_clair3:
	input:
		genome=config['genome_fasta'],
		genome_index=config['genome_fasta']+'.fai',
		indexed_bam=output_folder+'/sorted_bam_files/{sample}_output_sorted.bam.bai',
		sorted_bam=output_folder+'/sorted_bam_files/{sample}_output_sorted.bam',
	params:
		model_path=config['model_path'],
		snp_min_af=config['snp_min_af'],
		indel_min_af=config['indel_min_af'],
		clair3_folder=output_folder+'/clair3_samples/{sample}'
	output:
		clair3_gvcf=output_folder+'/clair3_samples/{sample}/merge_output.gvcf.gz'
	shell:
		'''
		run_clair3.sh \
		--bam_fn={input.sorted_bam} \
		--ref_fn={input.genome} \
		--threads=6 \
		--platform=ont \
		--model_path={params.model_path} \
		--output={params.clair3_folder} \
		--print_ref_calls \
		--include_all_ctgs \
		--gvcf \
		--no_phasing_for_fa \
		--sample_name={wildcards.sample} \
		--snp_min_af={params.snp_min_af} \
		--indel_min_af={params.indel_min_af}
		'''