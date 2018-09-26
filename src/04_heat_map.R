library(ggplot2)
#library(gplots)
library(data.table)
library(DESeq2)

dds_aru <- readRDS("output/dds_aru.Rds")
results_table <- readRDS("output/complete_results_table.Rds")
sample_info <- fread("data/sample_id_mapping.csv")

# transform read counts
vst <- varianceStabilizingTransformation(dds_aru, blind = FALSE)

# deseq pca
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(500)]
pca <- prcomp(t(assay(vst)[select, ]))
loadings <- pca$rotation

# get the highest contributing genes to PC1
loadings_dt <- data.table(loadings, keep.rownames = TRUE)
loadings_dt[, abs_pc1 := abs(PC1)]
setorder(loadings_dt, -abs_pc1)

number_to_plot <- 200
genes_to_plot <- loadings_dt[1:number_to_plot, unique(rn)]

# set plot data
lfc_pd <- results_table[denom == "Untreated" & gene_name %in% genes_to_plot]

# order the x axis
trt_order <- c("E. coli", "ARU", "Cyto")
lfc_pd[, numerator := factor(numerator, levels = trt_order)]

# order the y-axis
lfc_matrix <- as.matrix(data.frame(dcast(lfc_pd,
                 gene_name ~ numerator,
                 value.var = "log2FoldChange"),
           row.names = "gene_name"))
hc <- hclust(dist(lfc_matrix))
gene_order <- hc$labels[hc$order]
lfc_pd[, gene_name := factor(gene_name, levels = gene_order)]


# plot LFC for these genes
hs <- RColorBrewer::brewer.pal(6, "RdBu")
gp <- ggplot(lfc_pd, aes(x = numerator, y = gene_name, fill = log2FoldChange)) +
    theme_minimal(base_size = 6) +
    theme(axis.text.y = element_text(face = "italic", size = 1.5),
          legend.key.size = unit(3, "mm")) +
    scale_x_discrete(expand = c(0, 0)) +
    xlab(NULL) + ylab(NULL) +
    scale_fill_gradientn(colours = rev(hs),
                         guide = guide_colourbar(title = "L2FC")) +
    geom_raster()

ggsave("output/plots/heatmap.pdf", gp, width = 70, height = 110, units = "mm")


# 
# # generate vst data table
# vst_table_wide <- data.table(data.frame(assay(vst)), keep.rownames = TRUE)
# setnames(vst_table_wide, "rn", "gene_name")
# vst_table <- melt(vst_table_wide, id.vars = "gene_name",
#      variable.name = "sample_id",
#      value.name = "vst_counts")
# 
# # merge sample_info
# vst_with_info <- merge(vst_table, sample_info, by = "sample_id")
# 
# # sum
# vst_pd <- vst_with_info[gene_name %in% genes_to_plot,
#               .(mean_vst = mean(vst_counts)),
#               by = .(gene_name, treatment)]
# 
# # generate y-axis order
# vst_matrix <- as.matrix(data.frame(dcast(vst_pd, gene_name ~ treatment),
#                                    row.names = "gene_name"))
# hc <- hclust(dist(vst_matrix))
# gene_order <- hc$labels[hc$order]
# vst_pd[, gene_name := factor(gene_name, levels = gene_order)]
# 
# # specify x-axis order
# 
# # scale the vst_values
# scaled_vst <- t(scale(t(vst_matrix)))
# 
# # remove NAs
# keep <- apply(scaled_vst, 1, function(x) !any(is.na(x)))
# scaled_vst <- scaled_vst[keep,]
# 
# scaled_hc <- hclust(dist(scaled_vst, method = "minkowski"),
#                     method = "ward.D2")
# scaled_order <- scaled_hc$labels[scaled_hc$order]
# scaled_pd <- data.table(melt(scaled_vst))
# setnames(scaled_pd, c("Var1", "Var2", "value"),
#          c("gene_name", "treatment", "scaled_vst"))
# scaled_pd[, gene_name := factor(gene_name, levels = scaled_order)]
# 
# 
# 
# # draw with ggplot
# hs <- RColorBrewer::brewer.pal(6, "YlOrRd")
# 
# ggplot(vst_pd, aes(x = treatment, y = gene_name, fill = mean_vst)) +
#     theme_minimal() +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_fill_gradientn(colours = hs) +
#     xlab(NULL) + ylab(NULL) +
#     geom_raster()
# 
# ggplot(scaled_pd, aes(x = treatment, y = gene_name, fill = scaled_vst)) +
#     theme_minimal() +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_fill_gradientn(colours = hs) +
#     xlab(NULL) + ylab(NULL) +
#     geom_raster()
# 
# # draw with gplots
# heat_scale <- colorRampPalette(hs)
# 
# setcolorder(scaled_vst, c("E..coli", "X5.A.RU", "IL12.IL18", "Untreated"))
# 
# heatmap.2(scaled_vst,
#           Colv = c(2, 1, 3, 4),
#           dendrogram = "row",
#           scale = "none",
#           col = heat_scale,
#           trace = "none",
#           keysize = 1,
#           density.info = 'none',
#           margins = c(10, 5))
# 

