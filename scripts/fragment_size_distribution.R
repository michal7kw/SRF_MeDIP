library(Rsamtools)
library(ggplot2)

# Function to get fragment sizes from a BAM file
get_fragment_sizes <- function(bam_file) {
  param <- ScanBamParam(what = c("isize"))
  bam <- scanBam(bam_file, param = param)[[1]]
  abs(bam$isize[bam$isize != 0])  # Use absolute value and remove zero-length fragments
}

# Get list of BAM files
bam_files <- list.files(path = "path/to/bam/files", pattern = "*.bam", full.names = TRUE)

# Get fragment sizes for all BAM files
fragment_sizes <- lapply(bam_files, get_fragment_sizes)
names(fragment_sizes) <- basename(bam_files)

# Combine data for plotting
plot_data <- do.call(rbind, lapply(names(fragment_sizes), function(name) {
  data.frame(sample = name, size = fragment_sizes[[name]])
}))

# Create plot
p <- ggplot(plot_data, aes(x = size, color = sample)) +
  geom_density() +
  xlim(0, 1000) +
  theme_minimal() +
  labs(x = "Fragment Size", y = "Density", title = "Fragment Size Distribution")

# Save plot
ggsave("fragment_size_distribution.pdf", p, width = 10, height = 6)
