library(data.table)
library(VennDiagram)

results_table <- readRDS("output/complete_results_table.Rds")

# look for subsets of sig genes
results_table[contrast_name %in% c("ARU_Untreated", "Cyto_Untreated") &
                  padj < 0.05]

# extract all sig genes
de_by_contrast <- results_table[padj < 0.05,
                                .(de_genes = unique(gene_name)),
                                by = .(denom, numerator)]

# subset by contrast
unt_aru <- de_by_contrast[denom == "Untreated" & numerator == "ARU",
                          unique(de_genes)]
unt_cyt <- de_by_contrast[denom == "Untreated" & numerator == "Cyto",
                          unique(de_genes)]
unt_ec <- de_by_contrast[denom == "Untreated" & numerator == "E. coli",
                         unique(de_genes)]

# total DE genes
de_by_contrast[, length(unique(de_genes))]

# get subsets of the venn diagram, e.g.
unt_cyt[(!unt_cyt %in% unt_aru) &
            !(unt_cyt %in% unt_ec)]
vd_genes <- results_table[gene_name %in% c(unt_aru, unt_cyt, unt_ec)]
setorder(vd_genes, gene_name)

# draw venn diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list(
    "E. coli" = unt_ec,
    "5-A-RU" = unt_aru,
    "Cytokines" = unt_cyt),
    filename = NULL,
    fill = Set1,
    lty = "solid",
    lwd = 1,
    cex = 0.5,
    cat.cex = 0.5,
    fontfamily = 'Sans',
    cat.fontfamily = 'Sans',
    alpha = 0.5,
    margin = 0
)

# write the plot
cairo_pdf("output/plots/venn_diagram.pdf", width = 2.8, height = 2.4)
grid.newpage()
grid.draw(vd)
dev.off()
