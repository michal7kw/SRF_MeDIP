#!/usr/bin/env Rscript

# Load required libraries
library(GenomicRanges)
library(clusterProfiler)
library(org.Mm.eg.db)
library(KEGG.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(AnnotationDbi)

# Read configuration
config <- yaml::read_yaml("config.yaml")

# Read differential methylation results
diff_meth_file <- file.path(config$diff_meth_outdir, "differential_methylation_results.txt")
dmrs <- read.table(diff_meth_file, header = TRUE, stringsAsFactors = FALSE)

# Convert DMRs to GRanges object
dmr_ranges <- GRanges(
  seqnames = dmrs$chromosome,
  ranges = IRanges(start = dmrs$start, end = dmrs$end),
  strand = "*"
)

# Get nearest genes to DMRs
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)
nearest_genes <- unique(mcols(distanceToNearest(dmr_ranges, genes))$subject)

# Convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Mm.eg.db, 
                       keys = nearest_genes, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")

# Perform GO enrichment analysis
go_results <- enrichGO(gene = gene_symbols,
                       OrgDb = org.Mm.eg.db,
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)

# Perform KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(gene = gene_symbols,
                           organism = "mmu",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)

# Save results
write.csv(as.data.frame(go_results), 
          file = file.path(config$func_annot_outdir, "go_enrichment_results.csv"), 
          row.names = FALSE)

write.csv(as.data.frame(kegg_results), 
          file = file.path(config$func_annot_outdir, "kegg_enrichment_results.csv"), 
          row.names = FALSE)

# Generate plots
pdf(file = file.path(config$func_annot_outdir, "go_dotplot.pdf"), width = 10, height = 8)
dotplot(go_results, showCategory = 20)
dev.off()

pdf(file = file.path(config$func_annot_outdir, "kegg_dotplot.pdf"), width = 10, height = 8)
dotplot(kegg_results, showCategory = 20)
dev.off()

# Print summary
cat("Functional annotation analysis completed.\n")
cat("GO enrichment results saved to:", file.path(config$func_annot_outdir, "go_enrichment_results.csv"), "\n")
cat("KEGG enrichment results saved to:", file.path(config$func_annot_outdir, "kegg_enrichment_results.csv"), "\n")
cat("GO dotplot saved to:", file.path(config$func_annot_outdir, "go_dotplot.pdf"), "\n")
cat("KEGG dotplot saved to:", file.path(config$func_annot_outdir, "kegg_dotplot.pdf"), "\n")