library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Function to convert genomic coordinates to gene symbols
coords_to_genes <- function(df) {
  # Create GRanges object from DMR coordinates
  dmr_ranges <- GRanges(
    seqnames = df$chromosome,
    ranges = IRanges(start = df$start, end = df$end)
  )
  
  # Get gene coordinates from TxDb
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes <- genes(txdb)
  
  # Find overlaps between DMRs and genes
  overlaps <- findOverlaps(dmr_ranges, genes)
  
  # Get ENTREZ IDs for overlapping genes
  entrez_ids <- unique(names(genes)[subjectHits(overlaps)])
  
  # Convert ENTREZ IDs to gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db,
                        keys = entrez_ids,
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")
  
  return(gene_symbols)
}

# Read differential methylation results
tryCatch({
  diff_meth <- read.table(snakemake@input[[1]], header = TRUE, sep = "\t")
}, error = function(e) {
  stop("Error reading differential methylation results: ", e$message)
})

# Filter for significant DMRs (adjusted p-value < 0.05)
sig_dmrs <- diff_meth[diff_meth$adj_p_value < 0.05, ]

if(nrow(sig_dmrs) == 0) {
  warning("No significant DMRs found. Using relaxed threshold...")
  sig_dmrs <- diff_meth[diff_meth$p_value < 0.05, ]
}

# Get associated genes
genes <- coords_to_genes(sig_dmrs)

if(length(genes) == 0) {
  stop("No genes found associated with DMRs")
}

# Convert gene symbols to ENTREZ IDs for enrichment analysis
entrez_ids <- mapIds(org.Hs.eg.db,
                    keys = genes,
                    column = "ENTREZID",
                    keytype = "SYMBOL",
                    multiVals = "first")

entrez_ids <- unique(na.omit(entrez_ids))

# Perform GO enrichment analysis
go_enrich <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

# Perform KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Save results
tryCatch({
  # GO results
  write.csv(as.data.frame(go_enrich), 
            file = snakemake@output[["go_csv"]], 
            row.names = FALSE)
  
  # KEGG results
  write.csv(as.data.frame(kegg_enrich), 
            file = snakemake@output[["kegg_csv"]], 
            row.names = FALSE)
  
  # Create and save GO plot
  pdf(snakemake@output[["go_plot"]], width = 10, height = 8)
  print(dotplot(go_enrich, showCategory = 20) +
          ggtitle("GO Enrichment Analysis of DMRs") +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8)))
  dev.off()
  
  # Create and save KEGG plot
  pdf(snakemake@output[["kegg_plot"]], width = 10, height = 8)
  print(dotplot(kegg_enrich, showCategory = 20) +
          ggtitle("KEGG Pathway Enrichment Analysis of DMRs") +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 8)))
  dev.off()
}, error = function(e) {
  stop("Error saving results: ", e$message)
})

# Write summary statistics to log file
sink(snakemake@log[[1]])
cat("Functional Annotation Analysis Summary:\n")
cat("----------------------------------------\n")
cat("Total DMRs analyzed:", nrow(diff_meth), "\n")
cat("Significant DMRs:", nrow(sig_dmrs), "\n")
cat("Unique genes identified:", length(genes), "\n")
cat("GO terms enriched:", nrow(as.data.frame(go_enrich)), "\n")
cat("KEGG pathways enriched:", nrow(as.data.frame(kegg_enrich)), "\n")
sink()
