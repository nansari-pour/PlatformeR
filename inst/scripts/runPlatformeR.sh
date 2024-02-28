#!/bin/bash

# input CNA simulate file, set directories for chromosome fasta files, phased SNP files and output directory
cna_simulate_file="/full/path/to/cna_simulate_gw.txt"
out_dir="/full/path/to/out/directory/"
fasta_dir="/full/path/to/Ref_files/single_chr_hg38/"
phase_dir="/full/path/to/phase/input/files/"

# Tools
ART='/full/path/to/art_bin_MountRainier/art_illumina'
BWA='/package/bwa/0.7.17/bin/bwa'
SAM='/package/samtools/1.17/bin/samtools'

# Parameters
loh_haplotype="maternal"
gain_haplotype="maternal"
coverage=60
normal_coverage=30
ncores=8
simulated_purity=0.8

# R Module to load - provide R name/version according to your cluster system
module load R

# Get chromosome number from array job

# Check if we're in a job array environment
if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
  # SLURM system
  chrom="$SLURM_ARRAY_TASK_ID"
elif [[ -n "$SGE_TASK_ID" ]]; then
  # SGE system
  chrom="$SGE_TASK_ID"
else
  # Not in a known job array environment
  echo "Error: Not in a recognised job array environment (SLURM or SGE)."
  exit 1
fi

# echo chromosome
echo "PlatformeR execution for Chromosome $chrom"

# R script execution with all input arguments
Rscript /full/path/to/runPlatformeR.R $cna_simulate_file $chrom $out_dir $fasta_dir $phase_dir $ART $BWA $SAM $loh_haplotype $gain_haplotype $coverage $normal_coverage $ncores $simulated_purity
