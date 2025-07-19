run_clair3.sh \
--bam_fn=/home/alfred/other_people/bnsengim/np_amp_clair3_demo_07-07-25/scratch_space/PM2_output_sorted.bam \
--ref_fn=/home/alfred/other_people/bnsengim/np_amp_clair3_demo_07-07-25/Pf_3D7_genome/PlasmoDB-62_Pfalciparum3D7_Genome.fasta \
--threads=6 \
--platform=ont \
--model_path=/home/alfred/micromamba/envs/clair3/bin/models/r1041_e82_400bps_sup_v410 \
--output=finished_clair3 \
--print_ref_calls \
--include_all_ctgs \
--gvcf \
--no_phasing_for_fa \
--sample_name=PM2 \
--snp_min_af=0.001 \
--indel_min_af=0.001


#run_clair3.sh \
#--include_all_ctgs \
#--no_phasing_for_fa \
#--sample_name={wildcards.sample} \
#--gvcf \
#--snp_min_af=0.001 \
#--indel_min_af=0.001 \
#--print_ref_calls