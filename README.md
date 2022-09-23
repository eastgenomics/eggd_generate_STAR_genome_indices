# eggd_make_genome_indices

## What does this app do?
This app uses a .gtf file of transcripts and the Sentieon STAR tool to generate genome indices, which will annotate a b38 reference genome.

## What inputs are required for this app to run?
* `--sentieon_tar`: (file) Tarballed Sentieon package. Currently defaults to use Sentieon 202112.05
* `--reference_genome`:(file) Tarballed GRCh38 reference genome FASTA + index. Current defaults to use GRCh38.no_alt_analysis_set_chr_mask21.fasta-index.tar.gz in 001_Reference
* `--gtf_file`: (file) File providing gene transcript information. Currently defaults to the GENCODE gtf v41 (gencode.v41.annotation.gtf.gz)

## How does this app work?
eggd_make_genome_indices takes an input .gtf file of transcript data, and a reference genome. It creates genome indices.

## What does this app output?
eggd_make_genome_indices outputs a gzipped file of genome indices.

## Notes
* This app is not ready for production use

## This app was made by East GLH