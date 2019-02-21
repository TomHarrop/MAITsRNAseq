library(data.table)
library(DESeq2)

# load sample information
sample_info <- fread("data/sample_id_mapping.csv")

# load non-excel files
datafiles_txt <- list.files(path = "data",
                            pattern = "Chip.*.txt",
                            full.names = TRUE)
names(datafiles_txt) <- sub(".bcmatrix.txt", "", basename(datafiles_txt))
data_list <- lapply(datafiles_txt, fread)

# merge count data from 4 chips
MeltCountData <- function(count_data_table) {
    melt(count_data_table, id.vars = "Gene",
         measure.vars = grep("IonXpress_[[:digit:]]{3}",
                             names(count_data_table),
                             value = TRUE),
         variable.name = "sample_id",
         value.name = "count")
}
molten_data_list <- lapply(data_list, MeltCountData)
molten_data <- rbindlist(molten_data_list, idcol = "chip")

# merge sample information
data_with_info <- merge(molten_data,
                        sample_info,
                        by = "sample_id",
                        all.x = TRUE)
data_with_info[, donor := paste0("pt", donor)]
data_with_info[treatment == "5-A-RU", treatment := "ARU"]
data_with_info[treatment == "IL12+IL18", treatment := "Cyto"]

# generate a list of genes to keep
mean_cutoff <- 5
max_cutoff <- 10 # this is crude, but you can use whatever you like
keep_mean <- data_with_info[, mean(count), by = .(Gene, treatment)][
    V1 > mean_cutoff, unique(Gene)]
keep_genes <- keep_mean
# keep_max <- data_with_info[, max(count), by = Gene][
#     V1 > max_cutoff, unique(Gene)]
# keep_genes <- union(keep_mean, keep_max)

# generate sample table for DESeq2
GenerateSampleTable <- function(count_data_with_info) {
    sample_table <- unique(count_data_with_info,
                           by = c("sample_id", "chip", "donor", "treatment"))
    sample_table[, c("Gene", "count") := NULL]
    sample_data <- data.frame(sample_table, row.names = "sample_id")
    
    # convert columns to factors
    sample_data$donor <- factor(sample_data$donor)
    trt_levels <- c("Untreated", "E. coli", "ARU", "Cyto")
    sample_data$treatment <- factor(sample_data$treatment, levels = trt_levels)
    sample_data$chip <- factor(sample_data$chip)
    sample_data
}

# generate count data for DESeq2
GenerateCountTable <- function(count_data_with_info) {
    count_table <- dcast(count_data_with_info,
                         Gene ~ sample_id,
                         value.var = "count")
    # FIXME. Excel has converted gene names to dates. Remove them for now.
    count_table_filtered <- count_table[complete.cases(count_table)]
    count_data <- as.matrix(data.frame(count_table_filtered, row.names = "Gene"))
    count_data
}

# generating the DESeq2 object
GenerateDeseqObject <- function(count_data_with_info) {
    count_data <- GenerateCountTable(count_data_with_info)
    sample_data <- GenerateSampleTable(count_data_with_info)
    dds <- DESeqDataSetFromMatrix(countData = count_data,
                                  colData = sample_data,
                                  design = ~ treatment + donor)
    dds
}

# split into table containing ARU and table for other comparisons and generate
# deseq objects
dds <- GenerateDeseqObject(data_with_info)

# filter the dds
dds_filtered <- dds[keep_genes,]

# save DESeq2 object
saveRDS(dds, file = "output/dds.Rds")
saveRDS(dds_filtered, file = "output/dds_filtered.Rds")
