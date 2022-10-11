#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json

mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/reference_genome_fasta_and_index
mkdir /home/dnanexus/gtf

# Unpack tarred input files and decompress gzipped files
tar xvzf /home/dnanexus/in/sentieon_tar/sentieon-genomics-*.tar.gz -C /usr/local
tar xvzf /home/dnanexus/in/reference_genome_fasta_and_index/*tar.gz -C /home/dnanexus/reference_genome_fasta_and_index
gunzip /home/dnanexus/in/annotated_transcripts_gtf/*.gz

source /home/dnanexus/license_setup.sh # run license setup script

export SENTIEON_INSTALL_DIR=/usr/local/sentieon-genomics-*

SENTIEON_BIN_DIR=$(echo $SENTIEON_INSTALL_DIR/bin)

export PATH="$SENTIEON_BIN_DIR:$PATH"

# NUMBER_THREADS input to STAR needs the number of cores on the server node
# This can be extracted from the DNAnexus instance type
INSTANCE=$(dx describe --json $DX_JOB_ID | jq -r '.instanceType')  # Extract instance type

# Output file name will be formatted as follows: 
# ref_filename-gtf_filename-readlengthN
REFPATH=/home/dnanexus/in/reference_genome_fasta_and_index/
REFNAME=( "$REFPATH"*.fasta-index.tar.gz )
# If there is no file matching the pattern, REFNAME will = *.fasta-index.tar.gz
# Therefore can filter incorrect filenames based on the presence of an asterisk
if [[ "$REFNAME" == *"*"* ]];
    then { echo "Input reference file name not found; possibly missing .fasta-index.tar.gz suffix" >&2; exit 1; }
fi
GTFPATH=/home/dnanexus/in/annotated_transcripts_gtf/
GTFNAME=( "$GTFPATH"*.gtf )
# If there is no file matching the pattern, GTFNAME will = *.gtf
# Therefore can filter incorrect filenames based on the presence of an asterisk
if [[ "$GTFNAME" == *"*"* ]];
    then { echo "Input transcripts gtf file name not found; possibly missing .gtf suffix" >&2; exit 1; }
fi
CUT_REFERENCE=${REFNAME##*/}
CUT_REFERENCE=${CUT_REFERENCE%.fasta-index.tar.gz*}
CUT_GTF=${GTFNAME##*/}
CUT_GTF=${CUT_GTF%.gtf*}
let READLENGTH=${read_length_minus_one}+1
FILENAME=ref_${CUT_REFERENCE}-gtf_${CUT_GTF}-readlength${READLENGTH}

# Configure output directories with output filename 
mkdir /home/dnanexus/$FILENAME
mkdir -p /home/dnanexus/out/$FILENAME

## Generate genome indices
# Define input variables for STAR command
NUMBER_THREADS=${INSTANCE##*_x}
export REFERENCE=/home/dnanexus/reference_genome_fasta_and_index/*.fa  # Reference genome, standard GRCh38
GTF=/home/dnanexus/in/annotated_transcripts_gtf/*gtf  # Input .gtf annotation file
OUTPUT_DIR=/home/dnanexus/$FILENAME

# Run STAR command to generate genome indices
sentieon STAR --runThreadN ${NUMBER_THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${OUTPUT_DIR} \
    --genomeFastaFiles ${REFERENCE} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang ${read_length_minus_one}

# Tar and gzip output file
tar -czvf $FILENAME.tar.gz /home/dnanexus/$FILENAME
# Move to /out/ folder to allow output to be uploaded
mv $FILENAME.tar.gz /home/dnanexus/out/$FILENAME

output_indices=$(dx upload home/dnanexus/out/$FILENAME/$FILENAME.tar.gz --brief)
dx-jobutil-add-output genome_indices "$output_indices" --class=file
dx-upload-all-outputs
