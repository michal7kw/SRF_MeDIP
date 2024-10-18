library(MEDIPS)
library(GenomicRanges)

# Read in sample information
samples <- read.csv("sample_info.csv", stringsAsFactors = FALSE)

# Function to create MEDIPS SET object
create_medips_set <- function(bam_file, bs, chr.select) {
  MEDIPS.createSet(file = bam_file, BSgenome = bs, extend = 200, shift = 0, 
                   uniq = 1, window_size = 100, chr.select = chr.select)
}

# Create MEDIPS SET objects for all samples
medips_sets <- lapply(samples$bam_file, create_medips_set, 
                      bs = BSgenome.Hsapiens.UCSC.hg38, 
                      chr.select = paste0("chr", 1:22))

# Perform differential methylation analysis
diff_meth <- MEDIPS.meth(MSet1 = medips_sets[samples$condition == "treatment"], 
                         MSet2 = medips_sets[samples$condition == "control"], 
                         CSet = NULL, p.adj = "BH", diff.method = "edgeR", 
                         MeDIP = T, CNV = F, minRowSum = 10)

# Write results to file
write.table(diff_meth, file = "differential_methylation_results.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
