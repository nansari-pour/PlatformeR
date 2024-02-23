# PlatformeR
In silico whole-genome sequencing (WGS) alignment file (BAM) generation tool with simulated copy number aberrations (CNA)

# SLURM SYSTEM
sbatch --partition=YOUR_SYSTEM_PARTITION_NAME --ntasks=8 --mem=40G -a 1-22 ./PlatformeR.sh
