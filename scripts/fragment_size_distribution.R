# Add at the beginning of the script
log_file <- file(snakemake@log[[1]], open="wt")
sink(log_file)
sink(log_file, type="message")

# Load required libraries
library(ggplot2)
library(dplyr)
library(Rsamtools)

# Function to get fragment sizes from BAM file
get_fragment_sizes <- function(bam_file) {
    # Read BAM file
    param <- ScanBamParam(what=c("qwidth"))
    bam <- scanBam(bam_file, param=param)[[1]]
    
    # Get fragment sizes
    sizes <- bam$qwidth
    
    # Remove any NA values
    sizes <- sizes[!is.na(sizes)]
    
    # Return sizes between 0 and 1000 bp
    sizes <- sizes[sizes > 0 & sizes < 1000]
    
    return(sizes)
}

# Get input BAM files
bam_files <- snakemake@input

# Process each BAM file
fragment_sizes <- data.frame()
for(bam_file in bam_files) {
    tryCatch({
        sizes <- get_fragment_sizes(bam_file)
        if(length(sizes) > 0) {
            sample_name <- basename(bam_file)
            fragment_sizes <- rbind(fragment_sizes,
                                  data.frame(size=sizes,
                                           sample=rep(sample_name, length(sizes))))
        }
    }, error=function(e) {
        warning(sprintf("Error processing %s: %s", bam_file, e$message))
    })
}

# Create density plot if we have data
if(nrow(fragment_sizes) > 0) {
    p <- ggplot(fragment_sizes, aes(x=size, color=sample)) +
        geom_density(adjust=1.5) +
        theme_minimal() +
        labs(x="Fragment Size (bp)",
             y="Density",
             title="Fragment Size Distribution") +
        theme(legend.position="bottom",
              legend.text=element_text(size=6))
    
    # Save plot
    ggsave(snakemake@output[[1]], p, width=10, height=6)
} else {
    stop("No valid fragment size data found in any of the input BAM files")
}

# Close log file connections
sink(type="message")
sink()
