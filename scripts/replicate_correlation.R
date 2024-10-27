# Add at the beginning of the script
log_file <- file(snakemake@log[[1]], open="wt")
sink(log_file)
sink(log_file, type="message")

# Load required libraries
library(corrplot)
library(reshape2)
library(Rsamtools)
library(GenomicRanges)

# Function to get coverage in bins
get_coverage <- function(bam_file, bin_size=1000) {
    # Create bins across the genome
    bins <- tileGenome(seqinfo(BamFile(bam_file)), 
                      tilewidth=bin_size, 
                      cut.last.tile.in.chrom=TRUE)
    
    # Get coverage
    coverage <- countOverlaps(bins, 
                            BamFile(bam_file), 
                            param=ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE)))
    
    return(coverage)
}

# Get input BAM files
bam_files <- snakemake@input

# Process each BAM file
coverage_matrix <- matrix(nrow=0, ncol=length(bam_files))
colnames(coverage_matrix) <- basename(bam_files)

# Get coverage for each file
for(i in seq_along(bam_files)) {
    tryCatch({
        coverage <- get_coverage(bam_files[i])
        if(nrow(coverage_matrix) == 0) {
            coverage_matrix <- matrix(0, nrow=length(coverage), ncol=length(bam_files))
            colnames(coverage_matrix) <- basename(bam_files)
        }
        coverage_matrix[,i] <- coverage
    }, error=function(e) {
        warning(sprintf("Error processing %s: %s", bam_files[i], e$message))
    })
}

# Calculate correlation matrix if we have data
if(nrow(coverage_matrix) > 0) {
    cor_matrix <- cor(coverage_matrix, method="pearson")
    
    # Create correlation plot
    pdf(snakemake@output[[1]], width=10, height=10)
    corrplot(cor_matrix, 
            method="color", 
            type="upper", 
            tl.col="black", 
            tl.srt=45,
            addCoef.col="black")
    dev.off()
} else {
    stop("No valid coverage data found in any of the input BAM files")
}

# Close log file connections
sink(type="message")
sink()
