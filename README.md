# eggd_generate_STAR_genome_indices

## What does this app do?
This app uses a .gtf file of transcripts and the Sentieon STAR tool to generate genome indices, which will annotate a reference genome.

## What inputs are required for this app to run?
* `--sentieon_tar`: (file) Tarballed Sentieon package. Currently defaults to use Sentieon 202112.05
* `--reference_genome_fasta_and_index`:(file) Tarballed reference genome FASTA + index. Current defaults to use GRCh38.no_alt_analysis_set_chr_mask21.fasta-index.tar.gz in 001_Reference
* `--annotated_transcripts_gtf`: (file) File providing gene transcript information. Currently defaults to the GENCODE gtf v41 (gencode.v41.annotation.gtf.gz)
* `--read_length`: (int) The read length of data with which the genome indices will be used. Standard for Illumina instruments is 100; so the default is 100.

## How does this app work?
eggd_generate_STAR_genome_indices takes an input .gtf file of transcript data, and a reference genome. It uses Sentieon's STAR to create genome indices for use with STAR Fusion.

## What does this app output?
eggd_generate_STAR_genome_indices outputs a gzipped file of genome indices. The output file name is in the format `ref_{reference_genome_filename}-gtf_{gtf_filename}-readlength{read_length_number}`

## This app was made by East GLH