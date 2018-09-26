library(data.table)
library(ggplot2)

results_table <- readRDS("output/complete_results_table.Rds")

contrast_order <- c("E. coli_Untreated", "ARU_Untreated",
                    "Cyto_Untreated", "E. coli_Cyto", 
                    "ARU_Cyto", "E. coli_ARU")
results_table[, contrast_name := factor(contrast_name, levels = contrast_order)]

ggplot(results_table, aes(x = log2FoldChange,
                          y = -log10(padj),
                          colour = padj < 0.05)) +
    theme_grey(base_size = 16) +
    xlab("Log2-fold change") +
    ylab(expression(-log[10](italic(p)[adj]))) +
    scale_colour_manual(values = c("black", "red", NULL),
                        guide = FALSE) +
    xlim(c(-10, 10)) +
    facet_wrap(~ contrast_name) +
    geom_point(alpha = 0.5, size = 2)

