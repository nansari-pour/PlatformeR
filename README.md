# PlatformeR

An in silico whole-genome sequencing (WGS) copy number aberration (CNA) simulator in R

## Contents

1. [About](#about)
2. [Installation](#installation)
3. [Usage](#usage)
   - [Pipeline](#pipeline)
     - [SLURM system](#slurm-system)
     - [SGE system](#sge-system)
4. [Reference files](#reference-files)
5. [Input and Output](#input-and-output)
6. [FAQ](#faq)
7. [License](#license)

## About

This R package is an in silico genomics tool designed to generate WGS-based alignment files (BAM) with simulated CNA. This tool is intended to assist researchers in testing hypotheses and benchmarking CNA tools using simulated data.

PlatformeR is designed for efficient performance, typically completing per-chromosome simulations in less than 6 hours with 8 cores, ensuring fast generation of simulated WGS data.

## Installation

To install the package, you can use the `devtools` library in R (>3.0.2):

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
```
### Install your package from GitHub

```R
devtools::install_github("nansari-pour/PlatformeR")
```
### Dependencies

PlatformeR has been carefully designed so that it has no dependencies on other R packages.

However, it requires the following external tools for its functionality:

1. [art](https://www.niehs.nih.gov/research/resources/software/biostatistics/art) for generating FASTQ files based on the FASTA files generated by PlatformeR (version 2016.06.05 tested)

2. [bwa](https://github.com/lh3/bwa) for generating alignment BAM files from FASTQ files (v0.7.17 tested)

3. [samtools](http://www.htslib.org/) for handling SAM and BAM files (v1.17 tested)

## Usage

### Pipeline

The package includes Bash scripts which submit and run PlatformeR for simulating data and running jobs on both SLURM and SGE systems. The Bash (and respective R) script *templates* are located in the `inst/scripts` directory and must be updated based on the specific needs of the user.

#### SLURM system

To submit a job on a SLURM system, use the following command:

```bash
PlatformeR_submit=$(sbatch --partition=short --ntasks=8 --mem=40G --time=0-06:00:00 -a 1-22 ./runPlatformeR.sh | awk '{print $4}')
# Submit the genome-wide merging script with a dependency on the PlatformeR job
sbatch --partition=YOUR_SYSTEM_PARTITION_NAME --ntasks=8 --mem=20G --time=0-12:00:00 --dependency=afterok:${PlatformeR_submit} ./WGS_bam_generator.sh
```
Replace YOUR_SYSTEM_PARTITION_NAME with the appropriate SLURM partition.

#### SGE system

To submit a job on an SGE system, use the following command (**not tested**):

```bash
PlatformeR_submit=$(qsub -q YOUR_SYSTEM_QUEUE_NAME -pe shmem 8 -l h_rt=6:00:00,h_vmem=40G -t 1-22 ./runPlatformeR.sh | awk '{print $3}')
# Submit the genome-wide merging script with a dependency on the PlatformeR job
qsub -pe shmem 8 -l h_rt=12:00:00,h_vmem=20G -hold_jid ${PlatformeR_submit} ./WGS_bam_generator.sh
```
Replace YOUR_SYSTEM_QUEUE_NAME with the appropriate SGE queue.

**Note: Adjust the SLURM/SGE parameters based on your system's specifications but bear in mind PlatformeR will need ~40G of mem and around 5-6 hours to generate larger chromosomes.**

## Reference files

1. **Reference Genome Files:**
   - A set of chromosome-level hg38 reference genome files (FASTA format). These per-chromosome reference files are available at `http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/` which need to be indexed using **bwa**.

2. **Reference Phase Files:**
   - A set of reference phase files (TXT format) containing phased 1000G (hg38) heterozygous SNP allelic information per chromosome for individual NA12878. These files are available in `reference_files`.

## Input and Output

### Input:

1. **Input Tab-Delimited Text File:**
   - Include a tab-delimited text file specifying the parameters for the simulation of CNAs which includes five columns (with no headers) including 'chr', 'startpos', 'endpos', 'cna' (LOH or GAIN) and 'CCF' (0,1]. An example copy of this file is available in  `inst/extdata`. This file is unique to each simulation run.

Ensure that the input files are correctly formatted and that the necessary reference genome files are provided for accurate simulation.

### Output:

1. **WGS Tumour BAM File:**
   - Upon successful execution, PlatformeR will generate a Whole Genome Sequencing (WGS) tumor BAM file containing simulated reads based on the provided parameters.

2. **Optional: WGS Normal BAM File:**
   - If specified in the input parameters, an optional WGS normal BAM file with simulated reads may also be generated as the germline/normal control.

3. **Index Files:**
   - For each generated BAM file, corresponding index files (BAI) will be created to facilitate efficient data retrieval.

4. **Intermediary Files: For Advanced Users**
   - The working directory will contain intermediary files necessary for the simulation process. These include chromosome-level BAM files, stored in subdirectories named after chromosome identifiers (e.g., `chr1`, `chr2`, etc.). These intermediary files aid in the construction of the final BAM files and provide insights into the simulated data at various stages of the process.

The output files will be organized in the out directory (out_dir in the BASH scripts; working directory).

## FAQ

#### Q: Why is the tool called "PlatformeR"?

A: I have always been a big fan of platformer games and the more I analysed genome-wide copy number plots, the more they looked like the landscape of a platformer game, with the copy number aberrations resembling a platformed terrain. Hence, the name "PlatformeR" was chosen as a playful homage to platformer games.

#### Q: How can I contribute to PlatformeR development?

A: Contributions to PlatformeR are welcome! If you have ideas for new features, bug fixes, or improvements, please raise an issue or submit a pull request detailing the changes you've made and why they should be merged.
Your contributions will be reviewed, and feedback may be provided before merging. Thank you for helping improve PlatformeR!

## License

```
MIT License

Copyright (c) 2024 Naser Ansari-Pour

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
