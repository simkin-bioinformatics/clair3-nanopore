'''
parses the annotations of snpEff to output an AA table file for downstream use.
The current version does not handle targeted vs. nontargeted mutations. Basic
algorithm is as follows:
1. reads through a vcf file and stores sample columns
2. upon reaching a snpEff line, iterates first through all samples, then all 15
fields of snpEff. Because a single snpEff VCF line can contain multiple
mutations, any field that has remainder 0 in modulo 15 is treated as the
beginning of a new mutation. All fields are then stored in a nested dictionary,
keyed first by sample and then by mutation. format and clair3 output are also
captured
3. Search the gff file to look up gene names associated with gene IDs and use
snpEff annotations to fill in the 'header' portion of each mutation. Parse out
the clair3 'AD' field values into cov, ref, and alt, where ref is the first 'AD'
value, alt is the second, and coverage is the ref+alt (I don't use DP for
coverage because DP is always slightly bigger than summed AD values - probably
due to reads that clair3 considers 'sequencing error').
4. convert snpEff header fields and counts to AA tables.

The 6 AA table fields are:
Gene ID
Gene
Mutation Name
ExonicFunc
AA Change
Targeted

By comparison, the 15 fields snpEff outputs are:
'Allele 
 Annotation 
 Annotation_Impact 
 Gene_Name 
 Gene_ID 
 Feature_Type 
 Feature_ID 
 Transcript_BioType 
 Rank 
 HGVS.c 
 HGVS.p 
 cDNA.pos / cDNA.length 
 CDS.pos / CDS.length 
 AA.pos / AA.length 
 Distance 

The mapping of snpEff onto AA table is:
Gene_ID->Gene ID
Annotation->ExonicFunc,
gff lookup of Name associated with Gene_ID->Gene
HGVS.p->AA Change
Gene+AA Change->Mutation Name
Targeted - currently hardcoded as 'unspecified'
'''

#input_vcf='annotated_variants.vcf'
input_vcf='/home/alfred/other_people/bnsengim/np_amp_clair3_demo_07-07-25/clair3-nanopore/bienvenu_targeted_nanopore_alignment_no_filters/snp_Eff_output/annotated_variants.vcf'
input_gff='/home/alfred/other_people/bnsengim/np_amp_clair3_demo_07-07-25/snpEff/example_UCSC_gff_file/PlasmoDB-62_Pfalciparum3D7.gff'
targets_tsv='/home/alfred/other_people/bnsengim/np_amp_clair3_demo_07-07-25/targets.tsv'

def make_targets_dict(targets_tsv, gene_mappings):
	'''
	converts a targets.tsv file into a dictionary, such that the chromosome,
	reference allele, and alternate allele can be easily compared against
	entries from a vcf file. The bottom layer is a header that can be added to
	an output AA table. A later function will examine the VCF file lines and add
	sequencing depths.
	'''
	#[f'{ref+alt},{ref},{alt}', Gene_ID, gene_name, Mutation_Name, ExonicFunc, AA_change, targeted]
	targets_dict={}
	for line_number, line in enumerate(open(targets_tsv)):
		if line_number>0:
			line=line.strip().split('\t')
			chrom, pos, ID, ref, alt=line[0:5]
			mut=line[7]
			gene_ID=line[11]
			gene_name=gene_mappings[gene_ID]
			mutation_name=gene_name+'-'+mut
			targets_dict.setdefault(chrom, {})
			targets_dict[chrom].setdefault(pos, {})
			targets_dict[chrom][pos].setdefault(ref, {})
			targets_dict[chrom][pos][ref][alt]=[gene_ID, gene_name, mutation_name, 'missense_mutation', mut, 'Yes']
	return targets_dict

def extract_counts(labels, values):
	labels=labels.split(':')
	values=values.split(':')
	for label_number, label in enumerate(labels):
		if label=='AD':
			depths=values[label_number].split(',')
			break
	if len(depths)==2:
		ref, alt=depths
		ref=convert_count(ref)
		alt=convert_count(alt)
	elif len(depths)==1:
		ref=depths[0]
		ref=convert_count(ref)
		alt=0
	else:
		print('weird depths', labels, values)
	return ref, alt

def check_targeted(vcf_line, targets_dict, depth_dict, samples):
	'''
	checks a vcf_line to see if it's "targeted" using targets_dict and extracts
	reference and alternate counts for each sample if it is targeted. Outputs
	the resulting depths to depth_dict
	'''
	chrom, pos, ID, ref_allele, alt_allele=vcf_line[:5]
	mutation_number=1
	if chrom in targets_dict and pos in targets_dict[chrom] and ref_allele in targets_dict[chrom][pos]:
		for alt_allele in targets_dict[chrom][pos][ref_allele]:
			header_list=targets_dict[chrom][pos][ref_allele][alt_allele]
			for sample_number, sample in enumerate(samples):
				labels, sample_counts=vcf_line[8], vcf_line[9+sample_number]
				ref, alt=extract_counts(labels, sample_counts)
				depth_dict.setdefault(sample, {})
				depth_dict[sample][mutation_number]=[f'{ref+alt},{ref},{alt}']+header_list
			mutation_number+=1
	return depth_dict

def parse_vcf_file(input_vcf, targets_dict):
	'''
	searches vcf file for snpEff annotations and retrieves the mutation
	information and depths associated with the snpeff annotations
	'''
	ann_dict, target_depth_dict={},{}
	parsed_counter=0
	if input_vcf.endswith('.gz'):
		import gzip
		file_handle=gzip.open(input_vcf, mode='rt')
	else:
		file_handle=open(input_vcf)
	for line in file_handle:
		line=line.strip().split('\t')
		if line[0].startswith('#CHROM'):
			samples=line[9:]
		if len(line)>7 and line[7].split(';')[-1].startswith('ANN'):
			unparsed_snpeff=line[7].split('|')
			for sample_number, sample in enumerate(samples):
				for column_number, column in enumerate(unparsed_snpeff):
					parsed_column_number=column_number%15
					if parsed_column_number==0:
						parsed_counter+=1
						ann_dict.setdefault(sample, {})
						ann_dict[sample].setdefault(parsed_counter, [line[8]+';'+line[9+sample_number]])
					ann_dict[sample][parsed_counter].append(column)
		elif len(line)>7 and not line[0].startswith('#'):
			target_depth_dict=check_targeted(line, targets_dict, target_depth_dict, samples)
	return ann_dict, target_depth_dict

def convert_count(count):
	if count!='.':
		count=int(count)
	else:
		count=0
	return count

def parse_annotations(ann_dict, gene_mappings):
	'''
	extracts only columns of interest from the annotation dictionary
	'''
	protein_dict={}
	for sample in ann_dict:
		protein_dict.setdefault(sample, {})
		for mut in ann_dict[sample]:
			columns=ann_dict[sample][mut]
			if len(columns)==16 and columns[11]: #item 11 is only populated if the mutation is protein-coding
				vcf_fields=columns[0]
				Gene_ID=columns[5]
				gene_name=gene_mappings[columns[4]]
				AA_change=columns[11][2:]
				Mutation_Name=gene_name+'-'+AA_change
				ExonicFunc=columns[2]
				targeted='unspecified'
				labels, values=vcf_fields.split(';')
				ref, alt=extract_counts(labels, values)
				protein_dict[sample][mut]=[f'{ref+alt},{ref},{alt}', Gene_ID, gene_name, Mutation_Name, ExonicFunc, AA_change, targeted]
	return protein_dict


def grab_gene_mappings(input_gff):
	'''
	looks up gene names associated with gene IDs. Here is an example line from gff file to be parsed:
	ID=PF3D7_1343700;Name=Kelch13;description=kelch protein K13;ebi_biotype=protein_coding
	'''
	gene_mappings={}
	for line in open(input_gff):
		line=line.strip().split('\t')[-1].split(';')
		if line[0].startswith('ID=') and line[1].startswith('Name='):
			ID=line[0][3:]
			name=line[1][5:]
			gene_mappings[ID]=name
		elif line[0].startswith('ID='):
			ID=line[0][3:]
			name=ID
			gene_mappings[ID]=name
	return gene_mappings

def get_targeted_status(input_tsv, snpeff_dict):
	'''
	takes a targets.tsv file as input and looks up whether each mutation in a
	snpeff_dictionary is a 'targeted' mutation of interest or not.
	'''

def format_header(output_file, protein_dict):
	first_sample=list(protein_dict.keys())[0]
	header_names=['Gene ID', 'Gene', 'Mutation Name', 'ExonicFunc', 'AA Change', 'Targeted']
	for row_number, name in enumerate(header_names):
		output_line=[name]
		for mut in protein_dict[first_sample]:
			output_line.append(protein_dict[first_sample][mut][row_number+1])
		output_file.write(','.join(output_line)+'\n')

def output_tables(protein_dict, output_folder):
	for file_number, output_path in enumerate(['coverage_AA_table.csv', 'reference_AA_table.csv', 'alternate_AA_table.csv']):
		output_file=open(output_folder+'/'+output_path, 'w')
		format_header(output_file, protein_dict)
		for sample in protein_dict:
			output_line=[sample]
			for mut in protein_dict[sample]:
				output_line.append(protein_dict[sample][mut][0].split(',')[file_number])
			output_file.write(','.join(output_line)+'\n')

gene_mappings=grab_gene_mappings(input_gff)
targets_dict=make_targets_dict(targets_tsv, gene_mappings)
ann_dict, targeted_depth_dict=parse_vcf_file(input_vcf, targets_dict)
print('targeted dict is', targeted_depth_dict)
protein_dict=parse_annotations(ann_dict, gene_mappings)
output_tables(protein_dict, 'test_folder')