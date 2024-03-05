# PlatformeR
In silico whole-genome sequencing (WGS) alignment file (BAM) generation tool with simulated copy number aberrations (CNA)

## SLURM SYSTEM
sbatch --partition=YOUR_SYSTEM_PARTITION_NAME --ntasks=8 --mem=40G --time=0-06:00:00 -a 1-22 ./PlatformeR.sh

## SGE SYSTEM TEMPLATE (not tested)
qsub -q YOUR_SYSTEM_QUEUE_NAME -pe shmem 8 -l h_rt=6:00:00,h_vmem=40G -t 1-22 ./PlatformeR.sh

