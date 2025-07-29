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
		copied_snakefile=output_folder+'/snakemake_params/nanopore_variant_annotation.smk',
		snp_eff_vcf=output_folder+'/snp_Eff_output/annotated_variants.vcf'

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
		clair3_vcf=output_folder+'/clair3_samples/{sample}/merge_output.vcf.gz'
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

rule index_vcf:
	input:
		clair3_vcf=output_folder+'/clair3_samples/{sample}/merge_output.vcf.gz'
	output:
		vcf_index=output_folder+'/clair3_samples/{sample}/merge_output.vcf.gz.csi'
	shell:
		'bcftools index {input.clair3_vcf}'

rule merge_vcfs:
	input:
		clair3_vcfs=expand(output_folder+'/clair3_samples/{sample}/merge_output.vcf.gz', sample=sample_subset),
		vcf_indices=expand(output_folder+'/clair3_samples/{sample}/merge_output.vcf.gz.csi', sample=sample_subset)
	output:
		merged_vcf=output_folder+'/clair3_multisample/merged_multisample.vcf.gz'
	shell:
		'''
		bcftools merge \
			{input.clair3_vcfs} \
			--force-samples -O z -o {output.merged_vcf}
		'''

rule prepare_snp_eff_database:
	'''
	rearranges genome fasta and genome gff files for snp_eff. Specifically:
	1. makes a 'data' folder in the snpeff folder
	2. makes a subfolder in 'data' with your desired genome name
	3. move your gff and genome files to the newly created genome folder and
	renames your gff file as genes.gff and your genome as sequences.fa
	4. Adds a line to the config file underneath the data.dir line and puts the
	name of your genome folder plus .genome, followed by a description,
	formatted like this: your_desired_name.genome : your_description
	'''
	input:
		snp_eff_folder=config['snp_eff_folder'],
		genome_fasta=config['genome_fasta'],
		genome_gff=config['genome_gff']
	params:
		database=config['genome_database_name'],
		description=config['database_description']
	output:
		genome_fasta=config['snp_eff_folder']+'/data/'+config['genome_database_name']+'/sequences.fa',
		genome_gff=config['snp_eff_folder']+'/data/'+config['genome_database_name']+'/genes.gff',
		temp_config=output_folder+'/snp_eff_folder/snpEff.config'
	script:
		'scripts/prepare_snp_eff_database.py'

rule build_snp_eff:
	input:
		genome_fasta=config['snp_eff_folder']+'/data/'+config['genome_database_name']+'/sequences.fa',
	params:
		snp_eff_folder=config['snp_eff_folder'],
		database_name=config['genome_database_name']
	output:
		snp_eff_database=config['snp_eff_folder']+'/data/'+config['genome_database_name']+'/sequence.bin'
	shell:
		'java -Xmx8g -jar {params.snp_eff_folder}/snpEff.jar build -c {params.snp_eff_folder}/snpEff.config -noCheckCds -noCheckProtein {params.database_name}'

rule run_snp_eff:
	input:
		snp_eff_jar=config['snp_eff_folder']+'/snpEff.jar',
		merged_vcf=output_folder+'/clair3_multisample/merged_multisample.vcf.gz',
		snp_eff_database=config['snp_eff_folder']+'/data/'+config['genome_database_name']+'/sequence.bin'
	params:
		database_name=config['genome_database_name']
	output:
		snp_eff_stats=output_folder+'/snp_Eff_output/snpEff_summary.html',
		snp_eff_vcf=output_folder+'/snp_Eff_output/annotated_variants.vcf'
	shell:
		'''
		java -Xmx10g -jar {input.snp_eff_jar} {params.database_name} \
		{input.merged_vcf} \
		-stats {output.snp_eff_stats} >{output.snp_eff_vcf}
		'''

#rule parse_snp_eff: