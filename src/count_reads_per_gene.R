#!/usr/bin/env Rscript

library(BiocParallel)
library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)

# use 8 cores
register(MulticoreParam(workers = 8))

# read transcript db
gff_file <- "data/reference/GCF_000001405.37_GRCh38.p11_genomic.gff.gz"
gff <- makeTxDbFromGFF(gff_file, format = "gff3")
exons <- exonsBy(gff, "gene")

# read BAMFILE
bam <- BamFile("data/bam/IonXpress_017_R_2017_05_17_14_21_54_user_Liggins-Proton-104-JU_Transcriptome_Panel_CHIP_4_Auto_user_Liggins-Proton-104-JU_Transcriptome_Panel_CHIP_4_144.bam",
               index = "data/bam/IonXpress_017_R_2017_05_17_14_21_54_user_Liggins-Proton-104-JU_Transcriptome_Panel_CHIP_4_Auto_user_Liggins-Proton-104-JU_Transcriptome_Panel_CHIP_4_144.bam.bai",
               asMates = FALSE)

# count reads per gene
se <- summarizeOverlaps(features = exons,
                        reads = bam,
                        mode = "Union",
                        singleEnd = TRUE,
                        ignore.strand = TRUE)
