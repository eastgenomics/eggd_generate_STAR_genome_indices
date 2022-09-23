#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json

# resources folder will not exist on worker, so removed here.
mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/reference_genome
mkdir /home/dnanexus/output_indices
mkdir /home/dnanexus/gtf
mkdir -p /home/dnanexus/out/output_indices/

tar xvzf /home/dnanexus/in/sentieon_tar/sentieon-genomics-*.tar.gz -C /usr/local # unpack tar
tar xvzf /home/dnanexus/in/reference_genome/*tar.gz -C /home/dnanexus/reference_genome
gunzip /home/dnanexus/in/gtf_file/*.gz

source /home/dnanexus/license_setup.sh # run license setup script

export SENTIEON_INSTALL_DIR=/usr/local/sentieon-genomics-*

SENTIEON_BIN_DIR=$(echo $SENTIEON_INSTALL_DIR/bin)

export PATH="$SENTIEON_BIN_DIR:$PATH"

# NUMBER_THREADS input to STAR needs the number of cores on the server node
# This can be extracted from the DNAnexus instance type
INSTANCE=$(dx describe --json $DX_JOB_ID | jq -r '.instanceType')  # Extract instance type

## Generate genome indices
# Define input variables for STAR command
NUMBER_THREADS=${INSTANCE##*_x}
export REFERENCE=/home/dnanexus/reference_genome/*.fa  # Reference genome, standard GRCh38
GTF=/home/dnanexus/in/gtf_file/*gtf  # Input .gtf annotation file
READ_LENGTH_MINUS_1=99
OUTPUT_DIR=/home/dnanexus/output_indices

# Run STAR command to generate genome indices
sentieon STAR --runThreadN ${NUMBER_THREADS} \
    --runMode genomeGenerate \
    --genomeDir ${OUTPUT_DIR} \
    --genomeFastaFiles ${REFERENCE} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang ${READ_LENGTH_MINUS_1}

# Tar and gzip output file
tar -czvf output_indices.tar.gz /home/dnanexus/output_indices

# Move to /out/ folder to allow output to be uploaded
mv output_indices.tar.gz /home/dnanexus/out/output_indices

dx-upload-all-outputs
