library(DESeq2)
library(data.table)

# register parallel processors
BiocParallel::register(BiocParallel::MulticoreParam(8))

# load deseq objects
dds <- readRDS("output/dds_filtered.Rds")

# you can also filter on baseMean directly
# dds <- estimateSizeFactors(dds)
# norm_mean_by_gene <- rowMeans(counts(dds, normalized = TRUE))
# keep_genes <- names(norm_mean_by_gene[norm_mean_by_gene > 5])
# dds_filtered <- dds[keep_genes,]

# run deseq
dds <- DESeq(dds, parallel = TRUE)

# use a LRT
dds_lrt <- DESeq(dds, parallel = TRUE, test = "LRT", reduced = ~ donor)
res_lrt <- results(dds_lrt)
keep <- unique(rownames(subset(res_lrt, !is.na(padj))))
dds_filtered <- dds[keep, ]

# generate complete results table
ExtractResultsForContrast <- function(denominator_level,
                                      numerator_level) {
    my_deseq_object <- dds_filtered
    res <- results(my_deseq_object, contrast = c("treatment",
                                            numerator_level,
                                            denominator_level),
                   lfcThreshold = 1, alpha = 0.05)
    my_dt <- data.table(data.frame(res), keep.rownames = TRUE)
    setnames(my_dt, "rn", "gene_name")
    return(my_dt)
}

contrast_table <- data.table(
    denom = c("Untreated", "Untreated", "Untreated", "Cyto", "Cyto", "ARU"),
    numerator = c("E. coli", "ARU", "Cyto", "E. coli", "ARU", "E. coli")
)
results_table <- contrast_table[, ExtractResultsForContrast(
    denominator_level = denom,
    numerator_level = numerator),
    by = .(denom, numerator)]
results_table[, contrast_name := paste(numerator, denom, sep = "_")]

# setorder(results_table, baseMean)
# results_table[padj < 0.05]
# results_table[padj < 0.05 & baseMean < 5][, length(unique(gene_name))]
# results_table[padj < 0.05, length(unique(gene_name))]
# counts(dds, normalized = TRUE)["TM4SF4",]
# counts(dds, normalized = FALSE)["TM4SF4",]
# plotCounts(dds, gene = "TM4SF4", intgroup = "treatment")

# generate padj table
p_matrix <- dcast(results_table,
                  gene_name ~ contrast_name,
                  value.var = "padj")
lfc_matrix <- dcast(results_table,
                    gene_name ~ contrast_name,
                    value.var = "log2FoldChange")
fwrite(p_matrix, "output/p_matrix.csv")
fwrite(lfc_matrix, "output/lfc_matrix.csv")

# in case you want both tables together
merge(p_matrix,
      lfc_matrix,
      by = "gene_name",
      suffixes = c("_padj", "_log2FoldChange"))

# write the full results
fwrite(results_table, "output/complete_results_table.csv")
saveRDS(results_table, "output/complete_results_table.Rds")

# write results per contrast
lapply(results_table[, unique(contrast_name)], function(x)
    fwrite(results_table[contrast_name == x],
           paste0("output/contrasts/", x, ".csv")))
