#!/bin/bash

SCRIPT="Converting BCL to Fastq - v1.2 (Dec 23, 2021)"
AUTHOR="Copyright 2021 Usama Bakry (u.bakry@icloud.com)"

## Command line options
## -----------------------------------------------------------------------------
while getopts i:o:s:t:e: OPTION
do
case "${OPTION}"
in
# Input runfolder directory
i) INPUT=${OPTARG};;
# Parent output directory
o) OUTPUT=${OPTARG};;
# Sample sheet file
s) SAMPLESHEET=${OPTARG};;
# Threads
t) THREADS=${OPTARG};;
# Email
e) EMAIL=${OPTARG};;

esac
done
## -----------------------------------------------------------------------------

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] Input Run Directory  : ${INPUT}"
echo -e "[     INFO    ] Fastq Output Directory : ${OUTPUT}/$(basename "${INPUT}")"
echo -e "[     INFO    ] Sample Sheet : ${SAMPLESHEET}"
echo -e "[     INFO    ] Threads = ${THREADS}"
echo -e "[     INFO    ] Emails = ${EMAIL}"
## -----------------------------------------------------------------------------

## User confirmation
## -----------------------------------------------------------------------------
read -p "Continue (y/n)?" CHOICE
case "$CHOICE" in 
	y|Y ) 

## Main function
## -----------------------------------------------------------------------------                     
main () {

## Print pipeline info
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] ${SCRIPT}"
echo -e "[     INFO    ] ${AUTHOR}"
## -----------------------------------------------------------------------------  

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] Input Run Directory  : ${INPUT}"
echo -e "[     INFO    ] Fastq Output Directory : ${OUTPUT}/$(basename "${INPUT}")"
echo -e "[     INFO    ] Sample Sheet : ${SAMPLESHEET}"
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
OUTPUT="${OUTPUT}/$(basename "${INPUT}")/"

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

## Converting BCL to Fastq
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Converting BCL to Fastq..."

# bcl2fastq command
bcl2fastq -R $INPUT \
    -o $OUTPUT \
    --sample-sheet $SAMPLESHEET \
    --no-lane-splitting \
    -r $THREADS \
    -p $THREADS \
    -w $THREADS

echo -e "[      OK     ] Conversion is done.\n"

## -----------------------------------------------------------------------------

## Print end date/time
## -----------------------------------------------------------------------------                     
echo -e "[     END     ] $(date)\n"
## -----------------------------------------------------------------------------  

} ## End of main function
## -----------------------------------------------------------------------------                     

## Prepare output log file
## -----------------------------------------------------------------------------                     
LOG=$( { time main > "$(dirname "${OUTPUT}")"/bcl2fastq_output.log 2>&1; } 2>&1 )
echo -e "Duration:${LOG}" >> "$(dirname "${OUTPUT}")"/bcl2fastq_output.log 2>&1

mv "$(dirname "${OUTPUT}")"/bcl2fastq_output.log ${OUTPUT}/$(basename "${INPUT}")/
## -----------------------------------------------------------------------------                     

## Send email
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Sending email ..."

# ssmtp command
# ssmtp ${EMAIL} < ${OUTPUT}/bcl2fastq_output.log

# mail command
echo -e "[     INFO    ] ${SCRIPT}
[     INFO    ] ${AUTHOR}
                 
[     INFO    ] Input Run Directory  : ${INPUT}
[     INFO    ] Fastq Output Directory : ${OUTPUT}
[     INFO    ] Sample Sheet : ${SAMPLESHEET}
[     INFO    ] Threads = ${THREADS}
[     INFO    ] Emails = ${EMAIL}

Duration:${LOG}" | mail -s "$(basename "${INPUT}") - BCL2Fastq Conversion" -aFrom:Bioinformatics\ Operations\<bioinfo.ops@gmail.com\> -A ${OUTPUT}/$(basename "${INPUT}")/bcl2fastq_output.log ${EMAIL}

echo -e "[      OK     ] Email was sent.\n"

## -----------------------------------------------------------------------------

exit 0

;;

	n|N ) echo -e "[      OK     ] Process stopped.";;
	* ) echo -e   "[     ERROR   ] invalid";;
esac
