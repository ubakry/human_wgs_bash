#!/bin/bash

SCRIPT="Human Whole Genome Sequencing Data Analysis Pipeline (HWGS-PIPLINE) - v1.1 (Dec 26, 2021)"
AUTHOR="Copyright 2021 Usama Bakry (u.bakry@icloud.com)"

## Command line options
## -----------------------------------------------------------------------------
while getopts f:r:o:s:R:t:e: OPTION
do
case "${OPTION}"
in
# Forward fastq file
f) R1=${OPTARG};;
# Reverse fastq file
r) R2=${OPTARG};;
# Output directory
o) OUTPUT=${OPTARG};;
# sample name
s) SAMPLE_NAME=${OPTARG};;
# reference fasta file
R) REF=${OPTARG};;
# Threads
t) THREADS=${OPTARG};;
# Email
e) EMAIL=${OPTARG};;
esac
done
## -----------------------------------------------------------------------------

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] R1 Fastq File  : ${R1}"
echo -e "[     INFO    ] R2 Fastq File  : ${R2}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Threads = ${THREADS}"
echo -e "[     INFO    ] Emails = ${EMAIL}"
## -----------------------------------------------------------------------------

## User confirmation
## -----------------------------------------------------------------------------
read -p "Continue (y/n)?" CHOICE
case "$CHOICE" in 
	y|Y ) 

main () {

## Print pipeline info
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] ${SCRIPT}"
echo -e "[     INFO    ] ${AUTHOR}"
## -----------------------------------------------------------------------------  

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] R1 Fastq File  : ${R1}"
echo -e "[     INFO    ] R2 Fastq File  : ${R2}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Threads = ${THREADS}"
echo -e "[     INFO    ] Emails = ${EMAIL}\n"
## -----------------------------------------------------------------------------

## Print start date/time
## -----------------------------------------------------------------------------                     
echo -e "[    START    ] $(date)\n"
## -----------------------------------------------------------------------------  

## Create output directory
## -----------------------------------------------------------------------------                     
echo -e "[   PROCESS   ] Creating output directory..."

# Check if output directory is exist
if [ -d "${OUTPUT}" ] 
then
    echo -e "[    ERROR    ] The output directory already exists.\n"
    exit 0
else
    mkdir -p $OUTPUT/
fi

echo -e "[      OK     ] Output directory is ready on ${OUTPUT}\n"
## ----------------------------------------------------------------------------- 

## Quality control using fastp
## -----------------------------------------------------------------------------
time {
echo -e "[   PROCESS   ] Applying quality control..."

# Create qc subfolder
QC_DIR=${OUTPUT}/01_QC/
mkdir -p $QC_DIR

# Fastp command
fastp -i $R1 \
    -I $R2 \
    -o ${QC_DIR}/${SAMPLE_NAME}.R1.fq.gz \
    -O ${QC_DIR}/${SAMPLE_NAME}.R2.fq.gz \
    -j ${QC_DIR}/${SAMPLE_NAME}.json \
    -h ${QC_DIR}/${SAMPLE_NAME}.html \
    --thread $THREADS

echo -e "[      OK     ] QC is done.\n"
}
## -----------------------------------------------------------------------------

## Mapping against reference genome
## -----------------------------------------------------------------------------
time {
echo -e "[   PROCESS   ] Mapping against reference genome..."

# Create mapping subfolder
M_DIR=${OUTPUT}/02_Mapping/
mkdir -p $M_DIR

CHRS_DIR=${M_DIR}/Chrs/
mkdir -p $CHRS_DIR

# Index reference genome
#bwa index $REF

# BWA command
bwa mem $REF \
    ${QC_DIR}/${SAMPLE_NAME}.R1.fq.gz \
    ${QC_DIR}/${SAMPLE_NAME}.R2.fq.gz \
    -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}" \
    -t $THREADS > ${M_DIR}/${SAMPLE_NAME}.sam

# Convert and sort SAM file
samtools view -bS \
    ${M_DIR}/${SAMPLE_NAME}.sam \
    --threads $THREADS | \
samtools sort -o ${M_DIR}/${SAMPLE_NAME}.sorted.bam \
    --threads $THREADS

# Split BAM file and move files to Chrs directory
bamtools split -in ${M_DIR}/${SAMPLE_NAME}.sorted.bam -reference
mv ${M_DIR}/*.REF_chr*.bam ${CHRS_DIR}

echo -e "[      OK     ] Mapping is done.\n"
}
## -----------------------------------------------------------------------------

## Variant calling
## -----------------------------------------------------------------------------
time {
echo -e "[   PROCESS   ] Variant calling..."

# Zipping and indexing required files
# bgzip -c /SAN/DBs/UCSC_GRCh38/fasta/hg38_without_alt.fa > /SAN/DBs/UCSC_GRCh38/fasta/hg38_without_alt_bgzip.fa.gz
# gatk CreateSequenceDictionary -R /SAN/DBs/UCSC_GRCh38/fasta/hg38_without_alt_bgzip.fa.gz
# samtools faidx /SAN/DBs/UCSC_GRCh38/fasta/hg38_without_alt_bgzip.fa.gz

# Create variant calling subfolders
VC_DIR=${OUTPUT}/03_Variant_Calling/
mkdir -p $VC_DIR

GVCF_DIR=${VC_DIR}/GVCF/
mkdir -p $GVCF_DIR

HAPLO_CMD=""
CONVERT_CMD=""

BAMS=$(ls ${CHRS_DIR})
for BAM in $BAMS
do
    samtools index ${CHRS_DIR}/${BAM}

    BNAME_WITHOUT=$(basename "${BAM}" | cut -d. -f1,2,3)
    
    # GATK command
    HAPLO_CMD="${HAPLO_CMD} gatk --java-options \"-Xmx8g\" HaplotypeCaller -R /SAN/DBs/UCSC_GRCh38/fasta/hg38_without_alt_bgzip.fa.gz -I ${CHRS_DIR}/${BAM} -O ${VC_DIR}/${BNAME_WITHOUT}.g.vcf.gz -bamout ${VC_DIR}/${BNAME_WITHOUT}.bam -ERC GVCF &"
    CONVERT_CMD="${CONVERT_CMD} gatk --java-options \"-Xmx8g\" GenotypeGVCFs -R /SAN/DBs/UCSC_GRCh38/fasta/hg38_without_alt_bgzip.fa.gz -V ${VC_DIR}/${BNAME_WITHOUT}.g.vcf.gz -O ${VC_DIR}/${BNAME_WITHOUT}.vcf.gz &"
    
done

HAPLO_CMD="${HAPLO_CMD} wait && fg"
CONVERT_CMD="${CONVERT_CMD} wait && fg"

echo -e "${HAPLO_CMD}"
echo -e "${CONVERT_CMD}"

eval ${HAPLO_CMD}
eval ${CONVERT_CMD}

# Move GVCF files to GVCF folder
mv ${VC_DIR}/*.g.vcf.gz* ${GVCF_DIR}

# Concatenate VCF files to one VCF
bcftools concat --output ${VC_DIR}/${SAMPLE_NAME}.vcf.gz --output-type z --threads $THREADS $(ls ${VC_DIR}/*.vcf.gz)

echo -e "[      OK     ] Variant calling is done.\n"
}
## -----------------------------------------------------------------------------

## Variants annotation
## -----------------------------------------------------------------------------
time {
echo -e "[   PROCESS   ] Variants annotation..."

# Create annotation subfolder
VA_DIR=${OUTPUT}/04_Annotation/
mkdir -p $VA_DIR

# Unzip the vcf file
gunzip -k ${VC_DIR}/${SAMPLE_NAME}.vcf.gz

# VEP command
vep --input_file ${VC_DIR}/${SAMPLE_NAME}.vcf \
    --format vcf \
    --output_file $VA_DIR/${SAMPLE_NAME}.annot.vcf \
    --cache --warning_file $VA_DIR/${SAMPLE_NAME}.warning.log \
    --fork $THREADS \
    --everything \
    --verbose \
    --force_overwrite

echo -e "[      OK     ] Annotation is done.\n"
}
## -----------------------------------------------------------------------------

## Print end date/time
## -----------------------------------------------------------------------------                     
echo -e "[     END     ] $(date)\n"
## -----------------------------------------------------------------------------  

} ## End of main function
## -----------------------------------------------------------------------------                      

## Prepare output log file
## -----------------------------------------------------------------------------                     
LOG=$( { time main > "$(dirname "${OUTPUT}")"/${SAMPLE_NAME}.log 2>&1; } 2>&1 )
echo -e "Duration:${LOG}" >> "$(dirname "${OUTPUT}")"/${SAMPLE_NAME}.log 2>&1

# mv "$(dirname "${OUTPUT}")"/bcl2fastq_output.log ${OUTPUT}/$(basename "${INPUT}")/
## -----------------------------------------------------------------------------                     

## Send email
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Sending email ..."

# ssmtp command
# ssmtp ${EMAIL} < ${OUTPUT}/bcl2fastq_output.log

# mail command
echo -e "[     INFO    ] ${SCRIPT}
[     INFO    ] ${AUTHOR}
                 
[     INFO    ] R1 Fastq File  : ${R1}
[     INFO    ] R2 Fastq File  : ${R2}
[     INFO    ] Output Directory : ${OUTPUT}
[     INFO    ] Threads = ${THREADS}
[     INFO    ] Emails = ${EMAIL}

Duration:${LOG}" | mail -s "$(basename "${SAMPLE_NAME}") - WGS Data Analysis" -aFrom:Bioinformatics\ Operations\<bioinfo.ops@gmail.com\> -A $(dirname "${OUTPUT}")/${SAMPLE_NAME}.log ${EMAIL}

echo -e "[      OK     ] Email was sent.\n"

## -----------------------------------------------------------------------------


exit 0

;;

	n|N ) echo -e "[      OK     ] Process stopped.";;
	* ) echo -e   "[     ERROR   ] invalid";;
esac
