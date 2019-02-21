library(data.table)
library(VennDiagram)

results_table <- readRDS("output/complete_results_table.Rds")

# look for subsets of sig genes
results_table[contrast_name %in% c("ARU_Untreated", "Cyto_Untreated") &
                  padj < 0.05]
results_table[log2FoldChange > 0, direction := "up"]
results_table[log2FoldChange < 0, direction := "down"]

# extract all sig genes
de_by_contrast <- results_table[padj < 0.05,
                                .(de_genes = unique(gene_name)),
                                by = .(denom, numerator, direction)]

# subset by contrast
unt_aru_all <- de_by_contrast[denom == "Untreated" & numerator == "ARU",
                             unique(de_genes)]
unt_cyt_all <- de_by_contrast[denom == "Untreated" & numerator == "Cyto",
                             unique(de_genes)]
unt_ec_all <- de_by_contrast[denom == "Untreated" & numerator == "E. coli",
                            unique(de_genes)]

unt_aru_up <- de_by_contrast[denom == "Untreated" & numerator == "ARU" & direction == "up",
                          unique(de_genes)]
unt_cyt_up <- de_by_contrast[denom == "Untreated" & numerator == "Cyto" & direction == "up",
                          unique(de_genes)]
unt_ec_up <- de_by_contrast[denom == "Untreated" & numerator == "E. coli" & direction == "up",
                         unique(de_genes)]
unt_aru_down <- de_by_contrast[denom == "Untreated" & numerator == "ARU" & direction == "down",
                          unique(de_genes)]
unt_cyt_down <- de_by_contrast[denom == "Untreated" & numerator == "Cyto" & direction == "down",
                          unique(de_genes)]
unt_ec_down <- de_by_contrast[denom == "Untreated" & numerator == "E. coli" & direction == "down",
                         unique(de_genes)]

# why are there missing genes?
x <- intersect(unt_aru_all, unt_cyt_all)
y <- union(intersect(unt_aru_up, unt_cyt_up),
      intersect(unt_aru_down, unt_cyt_down))
wtf <- x[!x %in% y]
results_table[gene_name %in% wtf]

# total DE genes
de_by_contrast[, length(unique(de_genes))]

# get subsets of the venn diagram (NOT RUN)
# unt_cyt[(!unt_cyt %in% unt_aru) &
#             !(unt_cyt %in% unt_ec)]
# vd_genes <- results_table[gene_name %in% c(unt_aru, unt_cyt, unt_ec)]
# setorder(vd_genes, gene_name)

# draw venn diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")

# cex = number size

vd_all <- venn.diagram(x = list(
    "E. coli" = unt_ec_all,
    "5-A-RU" = unt_aru_all,
    "Cytokines" = unt_cyt_all),
    filename = NULL,
    fill = Set1,
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 1,
    cat.dist = 0.1,
    fontfamily = 'Sans',
    cat.fontfamily = 'Sans',
    alpha = 0.5,
    margin = 0.03)

vd_up <- venn.diagram(x = list(
    "E. coli" = unt_ec_up,
    "5-A-RU" = unt_aru_up,
    "Cytokines" = unt_cyt_up),
    filename = NULL,
    fill = Set1,
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 1,
    cat.dist = 0.1,
    fontfamily = 'Sans',
    cat.fontfamily = 'Sans',
    alpha = 0.5,
    margin = 0.03)

vd_down <- venn.diagram(x = list(
    "E. coli" = unt_ec_down,
    "5-A-RU" = unt_aru_down,
    "Cytokines" = unt_cyt_down),
    filename = NULL,
    fill = Set1,
    lty = "solid",
    lwd = 1,
    cex = 1,
    cat.cex = 1,
    cat.dist = 0.1,
    fontfamily = 'Sans',
    cat.fontfamily = 'Sans',
    alpha = 0.5,
    margin = 0.03)

# write the plot
cairo_pdf("output/plots/venn_diagram_up.pdf", width = 2.8, height = 2.4)
grid.newpage()
grid.draw(vd_up)
dev.off()

cairo_pdf("output/plots/venn_diagram_down.pdf", width = 2.3, height = 2.3)
grid.newpage()
grid.draw(vd_down)
dev.off()

cairo_pdf("output/plots/venn_diagram_all.pdf", width = 2.8, height = 2.4)
grid.newpage()
grid.draw(vd_all)
dev.off()
