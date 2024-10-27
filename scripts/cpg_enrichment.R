library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

# Function to calculate CpG enrichment
calculate_cpg_enrichment <- function(bam_file) {
  tryCatch({
    ms <- MEDIPS.createSet(file = bam_file, 
                          BSgenome = BSgenome.Hsapiens.UCSC.hg38, 
                          extend = 200, 
                          shift = 0, 
                          uniq = 1, 
                          window_size = 100, 
                          chr.select = paste0("chr", 1:22))
    
    enrich <- MEDIPS.CpGenrich(file = bam_file, 
                              BSgenome = BSgenome.Hsapiens.UCSC.hg38, 
                              extend = 200, 
                              shift = 0, 
                              uniq = 1)
    return(enrich)
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", bam_file, e$message))
    return(NA)
  })
}

# Get input BAM files from snakemake
bam_files <- snakemake@input[[1]]

# Calculate CpG enrichment for all BAM files
enrichment_results <- data.frame(
  file = bam_files,
  enrichment = sapply(bam_files, calculate_cpg_enrichment)
)

# Write results to file specified in snakemake output
write.table(enrichment_results, 
            file = snakemake@output[[1]], 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
