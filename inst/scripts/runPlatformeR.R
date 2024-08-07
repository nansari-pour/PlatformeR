# PlatformeR wrapper run

# add full Library path of PlatformeR installation to .libPaths()
.libPaths("/full/path/to/R/Library/where/PlatformeR/is/installed")

library(PlatformeR)

Args=commandArgs(trailingOnly = TRUE)
cna_simulate_file = toString(Args[1])
chrom = toString(Args[2])
chrom = ifelse(chrom == "X",chrom,as.integer(chrom))
print(paste("PlatformeR wrapper run for Chromosome",chrom))
out_dir = toString(Args[3])
fasta_dir = toString(Args[4])
phase_dir = toString(Args[5])
ART = toString(Args[6])
BWA = toString(Args[7])
SAM = toString(Args[8])
loh_haplotype = toString(Args[9])
gain_haplotype = toString(Args[10])
coverage = as.numeric(toString(Args[11]))
normal_coverage = as.numeric(toString(Args[12]))
ncores = as.integer(toString(Args[13]))
simulated_purity = as.numeric(toString(Args[14]))

print(paste("Coverage =",coverage,"Simulated Purity =",simulated_purity))

PlatformeR(cna_simulate_file = cna_simulate_file,
           chrom = chrom,
           out_dir = out_dir,
           fasta_dir = fasta_dir,
           phase_dir = phase_dir,
           art_bin = ART,
           coverage = coverage,
           read_length = 150,
           fragment_size = 400,
           fragment_size_sd = 10,
           tmp_dir = "TMP",
           bwa = BWA,
           ncores = ncores,
           samtools = SAM,
           loh_haplotype = loh_haplotype,
           gain_haplotype = gain_haplotype,
           simulated_purity = simulated_purity,
           chromX_haplotype = "paternal",
           generate_normal = TRUE,
           normal_coverage = normal_coverage,
           generate_diploid_normal = TRUE)
