library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

# Function to calculate CpG enrichment
calculate_cpg_enrichment <- function(bam_file) {
  ms <- MEDIPS.createSet(file = bam_file, BSgenome = BSgenome.Hsapiens.UCSC.hg38, 
                         extend = 200, shift = 0, uniq = 1, window_size = 100, 
                         chr.select = paste0("chr", 1:22))
  enrich <- MEDIPS.CpGenrich(file = bam_file, BSgenome = BSgenome.Hsapiens.UCSC.hg38, 
                             extend = 200, shift = 0, uniq = 1)
  return(enrich)
}

# Get list of BAM files
bam_files <- list.files(path = "path/to/bam/files", pattern = "*.bam", full.names = TRUE)

# Calculate CpG enrichment for all BAM files
enrichment <- sapply(bam_files, calculate_cpg_enrichment)

# Write results to file
write.table(enrichment, file = "cpg_enrichment.txt", sep = "\t", quote = FALSE)
