library(data.table)
library(ggplot2)

results_table <- readRDS("output/complete_results_table.Rds")

contrast_order <- c("E. coli_Untreated", "ARU_Untreated",
                    "Cyto_Untreated", "E. coli_Cyto", 
                    "ARU_Cyto", "E. coli_ARU")
results_table[, contrast_name := factor(contrast_name, levels = contrast_order)]

goi <- c("AREG", "VEGFA", "PDGFA")

# add labels and colours for GOI
results_table[gene_name %in% goi & padj < 0.05, gene_label := gene_name]
results_table[, pt_colour := 'a']
results_table[padj  < 0.05, pt_colour := 'b']
results_table[gene_name %in% goi & padj < 0.05, pt_colour := 'c']

gp <- ggplot(results_table,
             aes(x = log2FoldChange,
                          y = -log10(padj),
                          colour = pt_colour,
                          label = gene_label)) +
    theme_minimal(base_size = 8) +
    xlab("Log2-fold change") +
    ylab(expression(-log[10](italic(p)[adj]))) +
    scale_colour_manual(values = c("black", "red", "blue"),
                        guide = FALSE) +
    xlim(c(-10, 10)) +
    facet_wrap(~ contrast_name) +
    geom_point(shape = 16, alpha = 0.5, size = 2) +
    geom_text(colour = "blue", hjust = "outward", nudge_x = 0.25)


ggsave("output/plots/volcano_plot.pdf", gp, width = 10, height = 7.5, units = "in")
