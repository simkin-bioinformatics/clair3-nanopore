samtools faidx /home/alfred/other_people/bnsengim/np_amp_clair3_demo_07-07-25/Pf_3D7_genome/PlasmoDB-62_Pfalciparum3D7_Genome.fasta
samtools view -b -o PM2_output.bam PM2_output.sam
samtools sort -o PM2_output_sorted.bam PM2_output.bam
samtools index PM2_output_sorted.bam