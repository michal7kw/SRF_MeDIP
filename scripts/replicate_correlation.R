library(MEDIPS)
library(corrplot)
library(BSgenome.Hsapiens.UCSC.hg38)

# Function to create MEDIPS SET object
create_medips_set <- function(bam_file) {
  MEDIPS.createSet(file = bam_file, BSgenome = BSgenome.Hsapiens.UCSC.hg38, 
                   extend = 200, shift = 0, uniq = 1, window_size = 100, 
                   chr.select = paste0("chr", 1:22))
}

# Get list of BAM files
bam_files <- list.files(path = "path/to/bam/files", pattern = "*.bam", full.names = TRUE)

# Create MEDIPS SET objects for all samples
medips_sets <- lapply(bam_files, create_medips_set)

# Calculate correlation between samples
cor_matrix <- MEDIPS.correlation(MSets = medips_sets, plot = FALSE)

# Create correlation plot
pdf("replicate_correlation.pdf", width = 10, height = 10)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45)
dev.off()
