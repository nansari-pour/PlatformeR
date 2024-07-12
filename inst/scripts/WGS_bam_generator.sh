#!/bin/bash

# output directory
out_dir="/full/path/to/out/directory/"

# Tools
SAM='/package/samtools/1.17/bin/samtools'

# Parameters
chrom_names=$(echo $(seq 1 22) X)
echo $chrom_names
simulated_purity=0.8
ncores=8
tumourname="tumour"

# optional parameters - paired analysis
merge_normal="TRUE"
normalname="normal"

# R Module to load
module load R-cbrg/current

echo "PlatformeR genomewide/WGS BAM merge from individual chromosomes"

# R script execution with all input arguments
Rscript /full/path/to/WGS_bam_generator.R "$(echo $chrom_names)" $simulated_purity $out_dir $tumourname $SAM $ncores $merge_normal $normalname
