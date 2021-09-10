library(heatmaply)
library(tidyverse)

expression_heatmap <-
  function(values,
           gene_list,
           cols,
           row_font = 10,
           col_font = 10,
           cluster_rows = T,
           rnaseq_limits = NULL) {

    # Fixes the order in which I want the RNAseq to show
    # cols <- c("name",
    #           cols)


    # indices <- match(gene_list, rownames(values))
    # Skip missing
    # indices <- indices[!is.na(indices)]
    subgroup <- values[gene_list, cols]

    # rownames(subgroup) <- subgroup$name
    # subgroup$name <- NULL
    # colnames(subgroup) <- c("Ni_R1",
    #                         "Ni_R2",
    #                         "Ni_R3",
    #                         "Ni_EZH2i_R1",
    #                         "Ni_EZH2i_R2",
    #                         "Ni_EZH2i_R3",
    #                         "Pr_R1",
    #                         "Pr_R2",
    #                         "Pr_R3",
    #                         "Pr_EZH2i_R1",
    #                         "Pr_EZH2i_R2",
    #                         "Pr_EZH2i_R3")

    message("Missing : ", paste(values[!gene_list %in% rownames(subgroup), "gene"], sep = ","))

    # Make the RNA seq heatmap
    h <- heatmaply(log2(subgroup+1), Colv = F, Rowv = cluster_rows,
                   colorbar_xanchor='left', colorbar_yanchor='top', colorbar_xpos=1.1, colorbar_ypos=0.9, plot_method = "plotly",
                   fontsize_row = row_font, fontsize_col = col_font, limits = rnaseq_limits,
                   key.title = "log2(TPM + 1)")
    h
  }

top_upreg <- genes %>% filter(RNASeq_DS_EZH2i_vs_Ni_padj < 0.05 & RNASeq_DS_EZH2i_vs_Ni_log2FoldChange > 1) %>% arrange(desc(RNASeq_DS_EZH2i_vs_Ni_log2FoldChange))
top_upreg <- genes %>% filter(RNASeq_DS_EZH2i_vs_Ni_padj < 0.05 & RNASeq_DS_EZH2i_vs_Ni_log2FoldChange > 1) %>% arrange(RNASeq_DS_EZH2i_vs_Ni_padj)

top_upreg <- top_upreg[1:50, ]

datasets <- c("Kumar_2020", "Dong_2020")
datasets <- c("Collier_2017", "Kinoshita_2021", "Moody_2017", "Kumar_2020", "Collinson_2017", "Dong_2020", "Krendl_2017")
datasets <- c("Kumar_2020", "Krendl_2017")
counts_files <- file.path(params$rnaseqdir, datasets, "rsem.merged.gene_tpm.tsv")
counts_all <- lapply(counts_files, read_counts_file)

counts <- purrr::reduce(counts_all, full_join, by="gene_id")
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

cols_to_plot <- c("Kumar_2020_Naive_R1",
                  "Kumar_2020_Naive_R2",
                  "Kumar_2020_Naive_R3",
                  "Kumar_2020_Naive_EZH2i_R1",
                  "Kumar_2020_Naive_EZH2i_R2",
                  "Kumar_2020_Naive_EZH2i_R3",
                  "Kumar_2020_Primed_R1",
                  "Kumar_2020_Primed_R2",
                  "Kumar_2020_Primed_R3",
                  "Kumar_2020_Primed_EZH2i_R1",
                  "Kumar_2020_Primed_EZH2i_R2",
                  "Kumar_2020_Primed_EZH2i_R3",
                  "hESC_H9_undif_R1",
                  "hESC_H9_undif_R2",
                  "hESC_H9_8h_BMP4_R1",
                  "hESC_H9_8h_BMP4_R2",
                  "hESC_H9_16h_BMP4_R1",
                  "hESC_H9_16h_BMP4_R2",
                  "hESC_H9_24h_BMP4_R1",
                  "hESC_H9_24h_BMP4_R2",
                  "hESC_H9_48h_BMP4_R1",
                  "hESC_H9_48h_BMP4_R2",
                  "hESC_H9_72h_BMP4_R1",
                  "hESC_H9_72h_BMP4_R2"
                  )


genes_list <-
  c("EPAS1",
    "MSX2",
    "GATA3",
    "CLDN4",
    "GATA2",
    "IGF2",
    "CDX2",
    "SLC40A1",
    "KRT7",
    "FRZB",
    "CGA",
    "ERP27",
    "KRT23",
    "CGB5",
    "VGLL1",
    "ENPEP",
    "TP63"
  )

h <- expression_heatmap(counts, genes_list, cols_to_plot, row_font = 10, col_font = 10)

orca(h, "./output/fig5_lineage_markers_vs_Krendl_2017.svg", width = 9, height = 9, more_args = c('--disable-gpu'))

genes_list <- read.table("./data/other/Messmer_2019/Messmer_intermediate_up_top50.txt",
                         header = F)

h <- expression_heatmap(counts, genes_list$V1, cols_to_plot, row_font = 7, col_font = 7)

orca(h, "./output/fig4_messmer_top_50_up_vs_Krendl_2017.svg", width = 9, height = 9, more_args = c('--disable-gpu'))

genes_list <- read.table("./data/other/Messmer_2019/Messmer_intermediate_down_top50.txt",
                         header = F)

h <- expression_heatmap(counts, genes_list$V1, cols_to_plot, row_font = 7, col_font = 7)

orca(h, "./output/fig4_messmer_top_50_down_vs_Krendl_2017.svg", width = 9, height = 9, more_args = c('--disable-gpu'))

genes_list <- c("GATA2", "MSX2", "GATA3", "TFAP2A", "TFAP2C")

h <- expression_heatmap(counts, genes_list, cols_to_plot, row_font = 10, col_font = 10)

orca(h, "./output/fig5_Krendl_early_vs_Krendl_2017.svg", width = 9, height = 9, more_args = c('--disable-gpu'))



genes_list <- c("TBX3", "HAND1", "ANKRD1", "ISL1", "CDX2", "NR2F2", "DLX3", "ARID5B")

h <- expression_heatmap(counts, genes_list, cols_to_plot, row_font = 10, col_font = 10)

orca(h, "./output/fig5_Krendl_mid_vs_Krendl_2017.svg", width = 9, height = 9, more_args = c('--disable-gpu'))



genes_list <- c("EPAS1", "PPARG", "VGLL1", "MEIS1", "HOXB2", "BARX2", "TP63", "GCM1", "LCP1", "CEBPA", "BNC1", "GRHL1", "MEIS2", "TFAP2B", "SMARCA2")

h <- expression_heatmap(counts, genes_list, cols_to_plot, row_font = 10, col_font = 10)

orca(h, "./output/fig5_Krendl_late_vs_Krendl_2017.svg", width = 9, height = 9, more_args = c('--disable-gpu'))

