#!/bin/bash

### Converting BCL to Fastq
### Version 1.1 (Oct 17, 2021)
### Copyright 2021 Usama Bakry (u.bakry@icloud.com)

## Command line options
## -----------------------------------------------------------------------------
while getopts i:o:s:t: OPTION
do
case "${OPTION}"
in
# Input runfolder directory
i) INPUT=${OPTARG};;
# Output directory
o) OUTPUT=${OPTARG};;
# Sample sheet file
s) SAMPLESHEET=${OPTARG};;
# Threads
t) THREADS=${OPTARG};;
esac
done
## -----------------------------------------------------------------------------

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] Input Directory  : ${INPUT}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Sample Sheet : ${SAMPLESHEET}"
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
echo -e "[     INFO    ] Converting BCL to Fastq"
echo -e "[     INFO    ] Version 1.1 (Oct 17, 2021)"
echo -e "[     INFO    ] Copyright 2021 Usama Bakry (u.bakry@icloud.com)\n"
## -----------------------------------------------------------------------------  

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] Input Directory  : ${INPUT}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Sample Sheet : ${SAMPLESHEET}"
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

## Converting BCL to Fastq
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Converting BCL to Fastq..."

# bcl2fastq command
bcl2fastq -i $INPUT/Data/Intensities/BaseCalls \
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
echo -e "[     END     ] Starting date/time  : $(date)\n"
## -----------------------------------------------------------------------------  

exit 0
} > ${OUTPUT}/output.log 2>&1

;;

	n|N ) echo -e "[      OK     ] Process stopped.";;
	* ) echo -e   "[     ERROR   ] invalid";;
esac
