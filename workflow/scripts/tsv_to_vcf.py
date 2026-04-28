#!/usr/bin/env python3
import csv
import sys


def tsv_to_vcf(tsv_path, vcf_path):
    with open(tsv_path) as tsv_file, open(vcf_path, "w") as vcf_file:
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\n')
        vcf_file.write('##INFO=<ID=AA,Number=1,Type=String,Description="Amino acid change">\n')
        vcf_file.write('##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Gene ID">\n')
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        reader = csv.DictReader(tsv_file, delimiter="\t")
        for row in reader:
            info = f"GENE={row['gene_name']};AA={row['aminoacid_change']};GENE_ID={row['gene_id']}"
            vcf_file.write(
                f"{row['CHROM']}\t{row['POS']}\t.\t{row['REF']}\t{row['ALT']}\t.\t.\t{info}\n"
            )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.tsv> <output.vcf>", file=sys.stderr)
        sys.exit(1)
    tsv_to_vcf(sys.argv[1], sys.argv[2])
