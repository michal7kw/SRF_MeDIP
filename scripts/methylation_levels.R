library(MEDIPS)
library(GenomicRanges)
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

# Read in peak files
peak_files <- list.files(path = "path/to/peak/files", pattern = "*.narrowPeak", full.names = TRUE)
peaks <- lapply(peak_files, read.table, header = FALSE)
names(peaks) <- basename(peak_files)

# Combine all peaks
all_peaks <- do.call(rbind, lapply(names(peaks), function(name) {
  data.frame(chr = peaks[[name]]$V1, start = peaks[[name]]$V2, end = peaks[[name]]$V3, sample = name)
}))

# Create GRanges object from peaks
peak_ranges <- makeGRangesFromDataFrame(all_peaks)

# Calculate methylation levels in peak regions
meth_levels <- MEDIPS.meth(MSet1 = medips_sets, 
                           CSet = peak_ranges, 
                           CNV = FALSE, 
                           MeDIP = TRUE)

# Write results to file
write.table(meth_levels, file = "methylation_levels.txt", sep = "\t", quote = FALSE, row.names = FALSE)
