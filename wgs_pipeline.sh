#!/bin/bash

### Human Whole Genome Sequencing Data Analysis Pipeline (HWGS-PIPLINE)
### Version 1.0 (Sep 12, 2021)
### Copyright 2021 Usama Bakry and Ahmed ElHossieny

## Command line options
## -----------------------------------------------------------------------------
while getopts f:r:o:s:R:t: OPTION
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
esac
done
## -----------------------------------------------------------------------------

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] R1 Fastq File  : ${R1}"
echo -e "[     INFO    ] R2 Fastq File  : ${R2}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Threads = ${THREADS}"
## -----------------------------------------------------------------------------

## User confirmation
## -----------------------------------------------------------------------------
read -p "Continue (y/n)?" CHOICE
case "$CHOICE" in 
	y|Y ) 

time {

## Print pipeline info
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] Human Whole Genome Sequencing Data Analysis Pipeline (HWGS-PIPLINE)"
echo -e "[     INFO    ] Version 1.0 (Sep 12, 2021)"
echo -e "[     INFO    ] Copyright 2021 Usama Bakry and Ahmed ElHossieny\n"
## -----------------------------------------------------------------------------  

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] R1 Fastq File  : ${R1}"
echo -e "[     INFO    ] R2 Fastq File  : ${R2}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Threads = ${THREADS}\n"
## -----------------------------------------------------------------------------

## Print start date/time
## -----------------------------------------------------------------------------                     
echo -e "[    START    ] Starting date/time  : $(date)\n"
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

echo -e "[      OK     ] Mapping is done.\n"
}
## -----------------------------------------------------------------------------

## Variant calling
## -----------------------------------------------------------------------------
time {
echo -e "[   PROCESS   ] Variant calling..."

# Create variant calling subfolder
VC_DIR=${OUTPUT}/03_Variant_Calling/
mkdir -p $VC_DIR

# Zipping and indexing required files
bgzip -c /data1/ref/UCSC_GRCh38/hg38.fa > /data1/ref/UCSC_GRCh38/hg38_bgzip.fa.gz
gatk CreateSequenceDictionary -R /data1/ref/UCSC_GRCh38/hg38_bgzip.fa.gz
samtools faidx /data1/ref/UCSC_GRCh38/hg38_bgzip.fa.gz
samtools index ${M_DIR}/${SAMPLE_NAME}.sorted.bam

# GATK command
gatk --java-options "-Xmx400g" HaplotypeCaller  \
   -R /data1/ref/UCSC_GRCh38/hg38_bgzip.fa.gz \
   -I ${M_DIR}/${SAMPLE_NAME}.sorted.bam \
   -O ${VC_DIR}/${SAMPLE_NAME}.g.vcf.gz \
   -bamout ${VC_DIR}/${SAMPLE_NAME}.bam \
   -ERC GVCF
gatk --java-options "-Xmx400g" GenotypeGVCFs \
   -R /data1/ref/UCSC_GRCh38/hg38_bgzip.fa.gz \
   -V ${VC_DIR}/${SAMPLE_NAME}.g.vcf.gz \
   -O ${VC_DIR}/${SAMPLE_NAME}.vcf.gz

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
echo -e "[     END     ] Starting date/time  : $(date)\n"
## -----------------------------------------------------------------------------  

exit 0
} > output.log 2>&1

;;

	n|N ) echo -e "[      OK     ] Process stopped.";;
	* ) echo -e   "[     ERROR   ] invalid";;
esac
