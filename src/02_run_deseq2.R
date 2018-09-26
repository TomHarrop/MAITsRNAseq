library(DESeq2)
library(data.table)

# register parallel processors
BiocParallel::register(BiocParallel::MulticoreParam(8))

# load deseq objects
dds_aru <- readRDS("output/dds_aru.Rds")
dds_no_aru <- readRDS("output/dds_no_aru.Rds")

# run deseq
dds_aru <- DESeq(dds_aru, parallel = TRUE)
dds_no_aru <- DESeq(dds_no_aru, parallel = TRUE)

# generate complete results table
ExtractResultsForContrast <- function(denominator_level,
                                      numerator_level) {
    if (denominator_level == "ARU" || numerator_level == "ARU") {
        my_deseq_object <- dds_aru
    } else {
        my_deseq_object <- dds_no_aru
    }
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
