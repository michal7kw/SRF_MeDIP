library(MEDIPS)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(data.table)

# Function to create MEDIPS SET object
create_medips_set <- function(bam_file) {
  tryCatch({
    MEDIPS.createSet(
      file = bam_file,
      BSgenome = BSgenome.Hsapiens.UCSC.hg38,
      extend = 200,
      shift = 0,
      uniq = 1,
      window_size = 100,
      chr.select = paste0("chr", 1:22)
    )
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", bam_file, e$message))
    return(NULL)
  })
}

# Get input files from snakemake
bam_files <- snakemake@input[["bams"]]
peak_files <- snakemake@input[["peaks"]]

# Extract sample information
get_sample_info <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  list(
    condition = parts[2],  # DP, PP, or N
    type = parts[3],       # input or output
    replicate = gsub("r", "", parts[4])  # 1, 2, or 3
  )
}

# Organize samples
samples_df <- data.frame(
  file = bam_files,
  stringsAsFactors = FALSE
)

# Add sample information
sample_info <- lapply(samples_df$file, get_sample_info)
samples_df$condition <- sapply(sample_info, `[[`, "condition")
samples_df$type <- sapply(sample_info, `[[`, "type")
samples_df$replicate <- sapply(sample_info, `[[`, "replicate")

# Create MEDIPS SET objects
message("Creating MEDIPS SET objects...")
medips_sets <- lapply(bam_files, create_medips_set)
names(medips_sets) <- basename(bam_files)

# Remove any NULL entries from failed processing
medips_sets <- medips_sets[!sapply(medips_sets, is.null)]

# Read and combine peak files
message("Processing peak files...")
peaks <- lapply(peak_files, function(file) {
  tryCatch({
    peak_data <- fread(file, header = FALSE)
    parts <- strsplit(basename(file), "_")[[1]]
    condition <- parts[1]
    replicate <- gsub("_peaks.narrowPeak", "", parts[2])
    peak_data$condition <- condition
    peak_data$replicate <- replicate
    peak_data
  }, error = function(e) {
    warning(sprintf("Error reading peak file %s: %s", file, e$message))
    return(NULL)
  })
})

# Remove any NULL entries and combine
peaks <- do.call(rbind, peaks[!sapply(peaks, is.null)])

# Create GRanges object from peaks
peak_ranges <- GRanges(
  seqnames = peaks$V1,
  ranges = IRanges(start = peaks$V2, end = peaks$V3),
  condition = peaks$condition,
  replicate = peaks$replicate
)

# Calculate methylation levels for each condition
message("Calculating methylation levels...")
methylation_results <- list()

for (condition in unique(samples_df$condition)) {
  # Get treatment (output) samples for this condition
  treatment_files <- samples_df$file[samples_df$condition == condition & 
                                   samples_df$type == "output"]
  treatment_sets <- medips_sets[basename(treatment_files)]
  
  # Get control (input) samples for this condition
  control_files <- samples_df$file[samples_df$condition == condition & 
                                  samples_df$type == "input"]
  control_sets <- medips_sets[basename(control_files)]
  
  # Get peaks for this condition
  condition_peaks <- peak_ranges[peak_ranges$condition == condition]
  
  # Calculate methylation levels
  tryCatch({
    meth_levels <- MEDIPS.meth(
      MSet1 = treatment_sets,
      MSet2 = control_sets,
      CSet = condition_peaks,
      diff.method = "edgeR",
      MeDIP = TRUE,
      CNV = FALSE,
      minRowSum = 10
    )
    
    # Add condition information
    meth_levels$condition <- condition
    methylation_results[[condition]] <- as.data.frame(meth_levels)
  }, error = function(e) {
    warning(sprintf("Error calculating methylation levels for condition %s: %s", 
                   condition, e$message))
  })
}

# Combine results
message("Combining results...")
all_results <- do.call(rbind, methylation_results)

# Add additional information
all_results$peak_id <- paste(all_results$chr, all_results$start, all_results$end, sep = "_")

# Write results to file
message("Writing results to file...")
write.table(all_results,
            file = snakemake@output[[1]],
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Write summary to log file
sink(snakemake@log[[1]])
cat("Methylation Level Analysis Summary:\n")
cat("----------------------------------------\n")
cat("Total samples processed:", length(medips_sets), "\n")
cat("Total peaks analyzed:", length(peak_ranges), "\n")
cat("Conditions analyzed:", paste(unique(samples_df$condition), collapse = ", "), "\n")
cat("\nPer-condition summary:\n")
for (condition in unique(samples_df$condition)) {
  cat(sprintf("\n%s:\n", condition))
  cat(sprintf("  Peaks: %d\n", sum(peak_ranges$condition == condition)))
  cat(sprintf("  Samples: %d\n", sum(samples_df$condition == condition)))
}
sink()
