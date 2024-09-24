#!/bin/bash

# Tools
SAM='/package/samtools/1.17/bin/samtools'

# Parameters
chrom_names=$(seq 1 22) # autosomal chromosomes only for MSAI
echo $chrom_names
simulated_purity=0.8
ncores=8

# R Module to load - provide R name/version according to your cluster system
module load R

echo "PlatformeR genomewide/WGS BAM merge from individual chromosomes - TUMOUR 1 and PAIRED NORMAL"

# Merge chromosomes for TUMOUR1 and its paired NORMAL - output directory
out_dir="/full/path/to/out/directory/TUMOUR1/"
tumourname="tumour1"
# optional parameters - paired analysis
merge_normal="TRUE"
normalname="normal"

# R script execution with all input arguments
Rscript /full/path/to/WGS_bam_generator.R "$(echo $chrom_names)" $simulated_purity $out_dir $tumourname $SAM $ncores $merge_normal $normalname

echo "PlatformeR genomewide/WGS BAM merge from individual chromosomes - TUMOUR 2 ONLY"

# Merge chromosomes for TUMOUR2 ONLY - output directory
out_dir="/full/path/to/out/directory/TUMOUR2/"
tumourname="tumour2"
# optional parameters - paired analysis
merge_normal="FALSE"
normalname="none"

# R script execution with all input arguments
Rscript /full/path/to/WGS_bam_generator.R "$(echo $chrom_names)" $simulated_purity $out_dir $tumourname $SAM $ncores $merge_normal $normalname
