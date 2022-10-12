#!/bin/bash

set -exo pipefail #if any part goes wrong, job will fail

dx-download-all-inputs # download inputs from json

mkdir /home/dnanexus/genomeDir
mkdir /home/dnanexus/reference_genome_fasta_and_index

# Unpack tarred input files and decompress gzipped files
tar xvzf /home/dnanexus/in/sentieon_tar/sentieon-genomics-*.tar.gz -C /usr/local
tar xvzf /home/dnanexus/in/reference_genome_fasta_and_index/*tar.gz -C /home/dnanexus/reference_genome_fasta_and_index
gunzip /home/dnanexus/in/annotated_transcripts_gtf/*.gz

source /home/dnanexus/license_setup.sh # run license setup script

export SENTIEON_INSTALL_DIR=/usr/local/sentieon-genomics-*

SENTIEON_BIN_DIR=$(echo $SENTIEON_INSTALL_DIR/bin)

export PATH="$SENTIEON_BIN_DIR:$PATH"

# number_threads input to STAR needs the number of cores on the server node
# This can be extracted from the DNAnexus instance type
instance=$(dx describe --json $DX_JOB_ID | jq -r '.instanceType')  # Extract instance type

# Output file name will be formatted as follows: 
# ref_filename-gtf_filename-readlengthN
# First, extract the reference genome filename and the gtf filename from the path
refpath=/home/dnanexus/in/reference_genome_fasta_and_index/
refname=( "$refpath"*.fasta-index.tar.gz )
# If there is no file matching the pattern, refname will = *.fasta-index.tar.gz
# Therefore can filter incorrect filenames based on the presence of an asterisk
if [[ "$refname" == *"*"* ]];
    then { echo "Input reference file name not found; possibly missing .fasta-index.tar.gz suffix" >&2; exit 1; }
fi
gtfpath=/home/dnanexus/in/annotated_transcripts_gtf/
gtfname=( "$gtfpath"*.gtf )
# If there is no file matching the pattern, gtfname will = *.gtf
# Therefore can filter incorrect filenames based on the presence of an asterisk
if [[ "$gtfname" == *"*"* ]];
    then { echo "Input transcripts gtf file name not found; possibly missing .gtf suffix" >&2; exit 1; }
fi
cut_reference=${refname##*/}  # Cut after the final / in the path to get the filename
cut_reference=${cut_reference%.fasta-index.tar.gz*}  # Remove the suffix from the filename
# Repeat for GTF filename:
cut_gtf=${gtfname##*/} 
cut_gtf=${cut_gtf%.gtf*}
# Format output filename
filename=ref_${cut_reference}-gtf_${cut_gtf}-readlength${read_length}

# Create output directory with output filename 
output_dir=/home/dnanexus/$filename
mkdir $output_dir

## Generate genome indices
# Define input variables for STAR command
number_threads=${instance##*_x}
export REFERENCE=/home/dnanexus/reference_genome_fasta_and_index/*.fa  # Reference genome, standard GRCh38
gtf=/home/dnanexus/in/annotated_transcripts_gtf/*gtf  # Input .gtf annotation file
let read_length_minus_one=${read_length}-1

# Run STAR command to generate genome indices
sentieon STAR --runThreadN ${number_threads} \
    --runMode genomeGenerate \
    --genomeDir ${output_dir} \
    --genomeFastaFiles ${REFERENCE} \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang ${read_length_minus_one}

# Tar and gzip output file
tar -czvf $filename.tar.gz /home/dnanexus/$filename

# Upload output file
output_indices=$(dx upload /home/dnanexus/$filename.tar.gz --brief)
dx-jobutil-add-output genome_indices "$output_indices" --class=file