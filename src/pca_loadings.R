library(DESeq2)
library(data.table)
library(ggplot2)

dds <- readRDS("output/dds_filtered.Rds")
sample_info <- fread("data/sample_id_mapping.csv")


# transform read counts
vst <- varianceStabilizingTransformation(dds, blind = FALSE)

# deseq pca
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(500)]
pca <- prcomp(t(assay(vst)[select, ]))
loadings <- pca$rotation

# single-dimension pc plot
pc_long <- melt(data.table(pca$x, keep.rownames = TRUE), id.vars = "rn")
pc_long_pd <- merge(pc_long, sample_info, by.x = "rn", by.y = "sample_id")
ggplot(pc_long, aes(x = rn, y = value)) +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~variable) +
    geom_point()

ggplot(pc_long_pd, aes(x = treatment, y = value)) +
    facet_wrap(~variable) +
    geom_point()


# 2d PCA
pc_wide <- data.table(pca$x, keep.rownames = TRUE)
pc_wide_pd <- merge(pc_wide, sample_info, by.x = "rn", by.y = "sample_id")

ggplot(pc_wide_pd, aes(x = PC1, y = PC2, colour = treatment)) +
    geom_point()


# plot loadings
loading_pd <- melt(data.table(loadings, keep.rownames = TRUE),
                   id.vars = "rn",
                   variable.name = "component", value.name = "loading")

ggplot(loading_pd, aes(y = loading, x = rn)) +
    facet_wrap(~ component, scales = "free_y") +
    geom_point()

loading_pd <- data.table(loadings, keep.rownames = TRUE)

ggplot(loading_pd, aes(x = PC1, y = PC2, label = rn)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_text(alpha = 0.5)




