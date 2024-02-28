#!/bin/bash

# output directory
out_dir="/full/path/to/out/directory/"

# Tools
SAM='/package/samtools/1.17/bin/samtools'

# Parameters
chrom_names=$(seq 1 22)
echo $chrom_names
simulated_purity=0.8
ncores=8
samplename="tumour"

# optional parameters - paired analysis
merge_normal="TRUE"
normalname="normal"

# R Module to load
module load R-cbrg/current

echo "PlatformeR genomewide/WGS BAM merge from individual chromosomes"

# R script execution with all input arguments
Rscript /full/path/to/WGS_bam_generator.R "$(echo $chrom_names)" $simulated_purity $out_dir $samplename $SAM $ncores $merge_normal $normalname
