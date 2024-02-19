# Function to generate maternal and paternal chromosome fasta files per chromosome

prepare_chrom_fasta <- function(chrom,fasta_dir,phase_dir){
  
  print(chrom)
  
  # 1) read in het SNPs for the chromosome
  hetsnps=read.table(paste0(phase_dir,"NA12878_hg38_chr",chrom,"_snps.txt"),header = T, stringsAsFactors = F)
  
  # 2) read in the chromosome sequence
  fasta=readLines(paste0(fasta_dir,"chr",chrom,".fa"))
  length(fasta)
  chrom_fasta=paste(fasta[-1],collapse = "")
  chrom_fasta=toupper(chrom_fasta)
  names(chrom_fasta)=fasta[1]
  print(paste0("Chr",chrom," length= ",nchar(chrom_fasta)))
  chr_length <<- nchar(chrom_fasta)
  
  # 3) Generate paternal (ref) haplotype fasta file
  paternal_fasta_filename=paste0("chr",chrom,".paternal.fa")
  chrom_fasta_split=unlist(strsplit(chrom_fasta,split = ""))
  rm(chrom_fasta,fasta)
  chrom_fasta_split[hetsnps$POS] <- hetsnps$REF_ALLELE
  chrom_fasta_ref <<- paste0(chrom_fasta_split,collapse = "")
  writeLines(paste0(">chr",chrom,"-paternal"),paternal_fasta_filename)
  cat(chrom_fasta_ref, sep = "\n", file = paternal_fasta_filename, append = TRUE)
  print(paste0("Chr",chrom," paternal length = ",nchar(chrom_fasta_ref)))
  
  # 4) Generate maternal (alt) haplotype fasta file
  maternal_fasta_filename=paste0("chr",chrom,".maternal.fa")
  chrom_fasta_split[hetsnps$POS] <- hetsnps$ALT_ALLELE
  chrom_fasta_alt <<- paste0(chrom_fasta_split,collapse = "")
  rm(chrom_fasta_split)
  writeLines(paste0(">chr",chrom,"-maternal"),maternal_fasta_filename)
  cat(chrom_fasta_alt, sep = "\n", file = maternal_fasta_filename, append = TRUE)
  print(paste0("Chr",chrom," maternal length = ",nchar(chrom_fasta_alt)))
  
}
