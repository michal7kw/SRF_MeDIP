library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Read differential methylation results
diff_meth <- read.table("differential_methylation_results.txt", header = TRUE, sep = "\t")

# Extract significant differentially methylated regions (DMRs)
sig_dmrs <- diff_meth[diff_meth$adj.P.Value < 0.05, ]

# Convert genomic coordinates to gene IDs
genes <- unique(na.omit(select(org.Hs.eg.db, keys = sig_dmrs$gene_id, 
                               columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID))

# Perform GO enrichment analysis
go_enrich <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, ont = "BP", 
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Perform KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene = genes, organism = "hsa", 
                          pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(go_enrich), file = "go_enrichment_results.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_enrich), file = "kegg_enrichment_results.csv", row.names = FALSE)

# Create plots
pdf("go_dotplot.pdf", width = 10, height = 8)
dotplot(go_enrich, showCategory = 20)
dev.off()

pdf("kegg_dotplot.pdf", width = 10, height = 8)
dotplot(kegg_enrich, showCategory = 20)
dev.off()
