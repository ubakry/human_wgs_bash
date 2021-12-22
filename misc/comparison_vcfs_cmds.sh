# Unzip vcf file
gunzip -k <.vcf.gz>

# Rename the sample name in the vcf file
sed -i 's/<old_name>/<new_name>/g' <.vcf>

# Compress and indexing vcf file
bgzip -c <.vcf> > <.vcf.gz>
tabix -p vcf <.vcf.gz>

# Filter vcf file by chromosomes
bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrM,chrX,chrY <.vcf.gz> > <new.vcf>

# Get count of mutation in each chromosome 
cat <.vcf> | grep -v "^#" | cut -f1 | uniq -c

# Compare the two vcf files
time vcftools --vcf <ref_vcf> --diff <another_vcf> --out diff-site --diff-site
time vcftools --vcf <ref_vcf> --diff <another_vcf> --out diff-site-discordance --diff-site-discordance
time vcftools --vcf <ref_vcf> --diff <another_vcf> --out diff-discordance-matrix --diff-discordance-matrix

