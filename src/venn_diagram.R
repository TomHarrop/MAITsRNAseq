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

# get subsets of the venn diagram, e.g.
unt_cyt[(!unt_cyt %in% unt_aru) &
            !(unt_cyt %in% unt_ec)]
vd_genes <- results_table[gene_name %in% c(unt_aru, unt_cyt, unt_ec)]
setorder(vd_genes, gene_name)

# draw venn diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list(
    "5-A-RU" = unt_aru,
    "Cytokines" = unt_cyt,
    "E. coli" = unt_ec),
    filename = NULL,
    fill = Set1,
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 2,
    fontfamily = 'Sans',
    cat.fontfamily = 'Sans',
    alpha = 0.5,
    margin = 0
)
grid.newpage()
grid.draw(vd)
