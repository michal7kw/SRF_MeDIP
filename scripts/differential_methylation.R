library(MEDIPS)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Add at the beginning of the script
log_file <- file(snakemake@log[[1]], open="wt")
sink(log_file)
sink(log_file, type="message")

# Get input files from snakemake
bam_files <- snakemake@input[["bams"]]
peak_files <- snakemake@input[["peaks"]]

# Extract condition information from file names
get_condition <- function(filename) {
  # Extract condition from filename pattern "Medip_CONDITION_TYPE_rN"
  parts <- strsplit(basename(filename), "_")[[1]]
  return(parts[2])  # Returns DP, PP, or N
}

get_sample_type <- function(filename) {
  # Extract whether it's input or output
  parts <- strsplit(basename(filename), "_")[[1]]
  return(parts[3])  # Returns "input" or "output"
}

# Organize samples by condition
samples_df <- data.frame(
  file = bam_files,
  condition = sapply(bam_files, get_condition),
  type = sapply(bam_files, get_sample_type),
  stringsAsFactors = FALSE
)

# Function to create MEDIPS SET object
create_medips_set <- function(bam_file) {
  MEDIPS.createSet(
    file = bam_file,
    BSgenome = BSgenome.Hsapiens.UCSC.hg38,
    extend = 200,
    shift = 0,
    uniq = 1,
    window_size = 100,
    chr.select = paste0("chr", 1:22)
  )
}

# Create MEDIPS SET objects for all samples
medips_sets <- lapply(bam_files, function(bam) {
  tryCatch({
    create_medips_set(bam)
  }, error = function(e) {
    warning(sprintf("Error processing file %s: %s", bam, e$message))
    return(NULL)
  })
})

# Remove any NULL entries from failed processing
medips_sets <- medips_sets[!sapply(medips_sets, is.null)]

# Function to perform differential analysis between conditions
perform_diff_analysis <- function(condition1, condition2) {
  # Get treatment samples (output) for condition1
  treatment_files <- samples_df$file[samples_df$condition == condition1 & 
                                    samples_df$type == "output"]
  treatment_sets <- medips_sets[match(treatment_files, bam_files)]
  
  # Get control samples (output) for condition2
  control_files <- samples_df$file[samples_df$condition == condition2 & 
                                  samples_df$type == "output"]
  control_sets <- medips_sets[match(control_files, bam_files)]
  
  # Perform differential methylation analysis
  diff_meth <- MEDIPS.meth(
    MSet1 = treatment_sets,
    MSet2 = control_sets,
    CSet = NULL,
    p.adj = "BH",
    diff.method = "edgeR",
    MeDIP = TRUE,
    CNV = FALSE,
    minRowSum = 10
  )
  
  # Add comparison information
  diff_meth$comparison <- paste(condition1, "vs", condition2)
  return(diff_meth)
}

# Get unique conditions
conditions <- unique(samples_df$condition)

# Perform all pairwise comparisons
results_list <- list()
for (i in 1:(length(conditions)-1)) {
  for (j in (i+1):length(conditions)) {
    comparison_name <- paste(conditions[i], "vs", conditions[j])
    tryCatch({
      results_list[[comparison_name]] <- perform_diff_analysis(conditions[i], conditions[j])
    }, error = function(e) {
      warning(sprintf("Error in comparison %s: %s", comparison_name, e$message))
    })
  }
}

# Combine all results
all_results <- do.call(rbind, lapply(names(results_list), function(comparison) {
  res <- results_list[[comparison]]
  if (!is.null(res)) {
    data.frame(
      comparison = comparison,
      chromosome = res$chr,
      start = res$start,
      end = res$end,
      log2_fold_change = res$logFC,
      p_value = res$pvalue,
      adj_p_value = res$padj,
      stringsAsFactors = FALSE
    )
  }
}))

# Before writing results, add:
if (nrow(all_results) == 0) {
  stop("No differential methylation results found. Check input data and parameters.")
}

# Write results
write.table(all_results,
            file = snakemake@output[[1]],
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Add at the end of the script
sink(type="message")
sink()
