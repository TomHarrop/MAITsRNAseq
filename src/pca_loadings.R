library(DESeq2)
library(data.table)
library(ggplot2)

dds_aru <- readRDS("output/dds_aru.Rds")


# transform read counts
vst <- varianceStabilizingTransformation(dds_aru, blind = FALSE)

# deseq pca
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(500)]
pca <- prcomp(t(assay(vst)[select, ]))
loadings <- pca$rotation

# single-dimension pc plot
pc_long <- melt(data.table(pca$x, keep.rownames = TRUE), id.vars = "rn")
ggplot(pc_long, aes(x = rn, y = value)) +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~variable) +
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




