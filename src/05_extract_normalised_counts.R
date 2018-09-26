library(data.table)
library(DESeq2)

BiocParallel::register(BiocParallel::MulticoreParam(8))

# load  objects
dds_aru <- readRDS("output/dds_aru.Rds")
sample_info <- fread("data/sample_id_mapping.csv")
results_table <- fread("output/complete_results_table.csv")

# run deseq
dds_aru <- DESeq(dds_aru, parallel = TRUE)

# get count table
count_data_wide <- data.table(counts(dds_aru, normalized = TRUE),
                              keep.rownames = TRUE)
setnames(count_data_wide, "rn", "gene_id")
count_data_by_library <- melt(count_data_wide,
                              id.vars = "gene_id",
                              variable.name = "sample_id",
                              value.name = "normalised_count")

# make a table per library
count_data_long <- merge(count_data_by_library,
                         sample_info)

count_data <- dcast(count_data_long,
                    gene_id ~ treatment + donor,
                    value.var = "normalised_count")

# all significant genes
sig_genes_vs_untreated <- results_table[denom == "Untreated" &
                                            padj < 0.05, unique(gene_name)]

# write output
fwrite(count_data,
       "output/normalised_counts/all_genes.csv")
fwrite(count_data[gene_id %in% sig_genes_vs_untreated],
       "output/normalised_counts/sig_vs_untreated.csv")

lapply(results_table[, unique(contrast_name)],
       function(x)
       fwrite(count_data[gene_id %in%
                             results_table[contrast_name == x &
                                               padj < 0.05, unique(gene_name)]],
              paste0("output/normalised_counts/sig_", x, ".csv")))
