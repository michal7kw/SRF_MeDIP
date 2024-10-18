#!/usr/bin/env Rscript

# Load required libraries
required_packages <- c("GenomicRanges", "DESeq2", "rtracklayer", "yaml")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    cat(paste("Package", package, "is not installed. Please install it and try again.\n"))
    quit(status = 1)
  }
  library(package, character.only = TRUE)
}

# Error handling function
handle_error <- function(e) {
  cat("An error occurred:\n")
  cat(conditionMessage(e), "\n")
  cat("Call stack:\n")
  print(sys.calls())
  quit(status = 1)
}

tryCatch({
  cat("Starting script execution\n")
  
  # Read configuration
  cat("Reading configuration file\n")
  config_file <- "config.yaml"
  if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
  }
  
  # Read config file and ensure it ends with a newline
  config_lines <- readLines(config_file)
  config_text <- paste(c(config_lines, ""), collapse = "\n")
  
  config <- yaml::yaml.load(config_text)
  cat("Configuration loaded successfully\n")
  print(config)

  # Function to read MACS2 peak files
  read_peaks <- function(file_path) {
    cat("Reading file:", file_path, "\n")
    if (!file.exists(file_path)) {
      warning(paste("File does not exist:", file_path))
      return(NULL)
    }
    peaks <- tryCatch({
      import(file_path, format = "BED")
    }, error = function(e) {
      warning(paste("Error reading file:", file_path, "-", conditionMessage(e)))
      return(NULL)
    })
    if (is.null(peaks) || length(peaks) == 0) {
      warning(paste("File is empty or could not be read properly:", file_path))
      return(NULL)
    }
    cat("Successfully read", length(peaks), "peaks from", file_path, "\n")
    mcols(peaks)$score <- mcols(peaks)$signalValue
    return(peaks)
  }

  # Read all peak files
  cat("Reading peak files from directory:", config$peaks_outdir, "\n")
  peak_files <- list.files(config$peaks_outdir, pattern = "*_peaks.narrowPeak", full.names = TRUE)
  if (length(peak_files) == 0) {
    stop(paste("No peak files found in directory:", config$peaks_outdir))
  }
  cat("Found", length(peak_files), "peak files\n")
  
  peak_list <- lapply(peak_files, read_peaks)
  
  # Remove NULL entries (empty or unreadable files)
  peak_list <- peak_list[!sapply(peak_list, is.null)]
  cat("Successfully read", length(peak_list), "peak files\n")
  
  if (length(peak_list) == 0) {
    stop("No valid peak files found. Please check your input data.")
  }

  # Combine all peaks
  cat("Combining peaks\n")
  all_peaks <- unlist(GRangesList(peak_list))
  all_peaks <- reduce(all_peaks)
  cat("Combined", length(all_peaks), "unique peak regions\n")

  # Create count matrix
  cat("Creating count matrix\n")
  sample_names <- gsub("_peaks.narrowPeak", "", basename(peak_files))
  count_matrix <- matrix(0, nrow = length(all_peaks), ncol = length(sample_names))
  colnames(count_matrix) <- sample_names

  for (i in 1:length(peak_list)) {
    cat("Processing sample", i, "of", length(peak_list), "\n")
    overlaps <- findOverlaps(all_peaks, peak_list[[i]])
    count_matrix[queryHits(overlaps), i] <- mcols(peak_list[[i]])$score[subjectHits(overlaps)]
  }

  # Create sample information
  cat("Creating sample information\n")
  sample_info <- data.frame(
    condition = gsub("_.*", "", sample_names),
    replicate = gsub(".*_r", "", sample_names)
  )
  print(sample_info)

  # Check for sufficient number of conditions
  unique_conditions <- unique(sample_info$condition)
  if (length(unique_conditions) < 2) {
    stop("At least two different conditions are required for differential analysis.")
  }

  # Create DESeq2 object
  cat("Creating DESeq2 object\n")
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = sample_info,
                                design = ~ condition)

  # Run DESeq2
  cat("Running DESeq2\n")
  dds <- DESeq(dds)

  # Perform pairwise comparisons
  cat("Performing pairwise comparisons\n")
  comparisons <- combn(unique_conditions, 2, simplify = FALSE)

  results_list <- lapply(comparisons, function(comp) {
    cat("Comparing", comp[1], "vs", comp[2], "\n")
    res <- results(dds, contrast = c("condition", comp[1], comp[2]))
    res <- res[order(res$padj), ]
    return(res)
  })

  # Ensure output directory exists
  dir.create(config$diff_meth_outdir, showWarnings = FALSE, recursive = TRUE)

  # Write results
  cat("Writing results\n")
  write_results <- function(res, comparison, output_dir) {
    filename <- file.path(output_dir, paste0("diff_meth_", comparison[1], "_vs_", comparison[2], ".csv"))
    write.csv(as.data.frame(res), file = filename)
    cat("Results written to", filename, "\n")
  }

  mapply(write_results, results_list, comparisons, MoreArgs = list(output_dir = config$diff_meth_outdir))

  # Combine all results
  cat("Combining all results\n")
  all_results <- do.call(rbind, lapply(seq_along(results_list), function(i) {
    res <- as.data.frame(results_list[[i]])
    res$comparison <- paste(comparisons[[i]][1], "vs", comparisons[[i]][2])
    res$region <- rownames(res)
    return(res)
  }))

  # Write combined results
  cat("Writing combined results\n")
  write.table(all_results, 
              file = file.path(config$diff_meth_outdir, "differential_methylation_results.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Generate summary
  cat("Generating summary\n")
  summary_data <- lapply(seq_along(results_list), function(i) {
    res <- results_list[[i]]
    comp <- comparisons[[i]]
    data.frame(
      comparison = paste(comp[1], "vs", comp[2]),
      total_regions = nrow(res),
      significant_regions = sum(res$padj < 0.05, na.rm = TRUE),
      upregulated = sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
      downregulated = sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
    )
  })

  summary_table <- do.call(rbind, summary_data)
  write.csv(summary_table, file = file.path(config$diff_meth_outdir, "differential_methylation_summary.csv"), row.names = FALSE)

  # Print summary
  cat("Differential methylation analysis completed.\n")
  cat("Results saved to:", file.path(config$diff_meth_outdir, "differential_methylation_results.txt"), "\n")
  cat("Summary saved to:", file.path(config$diff_meth_outdir, "differential_methylation_summary.csv"), "\n")
  print(summary_table)

}, error = handle_error)