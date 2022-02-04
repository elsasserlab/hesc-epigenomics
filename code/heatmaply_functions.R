library(heatmaply)
library(tidyverse)

combined_heatmap <-
  function(values,
           gene_list,
           row_font = 10,
           col_font = 10,
           cluster_rows = T,
           rnaseq_limits = NULL,
           k27m3_limits = NULL,
           k4m3_limits = NULL,
           ub_limits = NULL) {

    # Fixes the order in which I want the RNAseq to show
    cols <- c("name",
              "RNASeq_TPM_Ni_R1",
              "RNASeq_TPM_Ni_R2",
              "RNASeq_TPM_Ni_R3",
              "RNASeq_TPM_Ni_EZH2i_R1",
              "RNASeq_TPM_Ni_EZH2i_R2",
              "RNASeq_TPM_Ni_EZH2i_R3",
              "RNASeq_TPM_Pr_R1",
              "RNASeq_TPM_Pr_R2",
              "RNASeq_TPM_Pr_R3",
              "RNASeq_TPM_Pr_EZH2i_R1",
              "RNASeq_TPM_Pr_EZH2i_R2",
              "RNASeq_TPM_Pr_EZH2i_R3")

    indices <- match(gene_list, values$name)
    # Skip missing
    indices <- indices[!is.na(indices)]
    subgroup <- values[indices, cols]

    rownames(subgroup) <- subgroup$name
    subgroup$name <- NULL
    colnames(subgroup) <- c("Ni_R1",
                            "Ni_R2",
                            "Ni_R3",
                            "Ni_EZH2i_R1",
                            "Ni_EZH2i_R2",
                            "Ni_EZH2i_R3",
                            "Pr_R1",
                            "Pr_R2",
                            "Pr_R3",
                            "Pr_EZH2i_R1",
                            "Pr_EZH2i_R2",
                            "Pr_EZH2i_R3")

    message("Missing : ", paste(setdiff(gene_list, rownames(subgroup)), sep = ",", collapse = ","))
    # Make the RNA seq heatmap
    h <- heatmaply(log2(subgroup+1), Colv = F, Rowv = cluster_rows,
                   colorbar_xanchor='left', colorbar_yanchor='top', colorbar_xpos=1.1, colorbar_ypos=0.9, plot_method = "plotly",
                   fontsize_row = row_font, fontsize_col = col_font, limits = rnaseq_limits,
                   key.title = "log2(TPM + 1)")
    # Get the order
    genes_ordered <- h$x$layout$yaxis$ticktext

    # Make the H3k27m3 heatmap with that order
    k27_cols <- c("name",
                  "H3K27m3_Ni_mean_cov",
                  "H3K27m3_Ni_EZH2i_mean_cov",
                  "H3K27m3_Pr_mean_cov",
                  "H3K27m3_Pr_EZH2i_mean_cov")

    k27_subgroup <- values[indices, k27_cols]
    rownames(k27_subgroup) <- k27_subgroup$name
    k27_subgroup$name <- NULL
    colnames(k27_subgroup) <- c("H3K27m3_Ni" ,
                                "H3K27m3_Ni_EZH2i",
                                "H3K27m3_Pr",
                                "H3K27m3_Pr_EZH2i")

    k27_colors <- colorRampPalette(c("white", "#3e5aa8"))(n = 1000)

    k27 <- heatmaply(k27_subgroup[rev(genes_ordered), ], Rowv = F, Colv = F,
                     colors = k27_colors,
                     fontsize_row = row_font, fontsize_col = col_font,
                     colorbar_xanchor='left',colorbar_yanchor='top', colorbar_xpos=1.1, colorbar_ypos=0.6, plot_method = "plotly",
                     key.title = "H3K27m3", limits = k27m3_limits)

    k4_cols <- c("name",
                 "H3K4m3_Ni_mean_cov",
                 "H3K4m3_Ni_EZH2i_mean_cov",
                 "H3K4m3_Pr_mean_cov",
                 "H3K4m3_Pr_EZH2i_mean_cov")

    k4_colors <- colorRampPalette(c("white", "#b64c28"))(n = 1000)

    k4_subgroup <- values[indices, k4_cols]

    rownames(k4_subgroup) <- k4_subgroup$name
    k4_subgroup$name <- NULL
    colnames(k4_subgroup) <- c("H3k4m3_Ni" ,
                               "H3k4m3_Ni_EZH2i",
                               "H3k4m3_Pr",
                               "H3k4m3_Pr_EZH2i")

    k4_labels <- k4_subgroup[rev(genes_ordered), ]

    k4 <- heatmaply(k4_subgroup[rev(genes_ordered), ], Rowv = F, Colv = F,
                    colors = k4_colors,
                    fontsize_row = row_font, fontsize_col = col_font,
                    colorbar_xanchor='left',colorbar_yanchor='top', colorbar_xpos=1.1, colorbar_ypos=0.3, colorbar_thickness = 30, plot_method = "plotly",
                    key.title = "H3K4m3", limits = k4m3_limits)


    ub_cols <- c("name",
                 "H2Aub_Ni_mean_cov",
                 "H2Aub_Ni_EZH2i_mean_cov",
                 "H2Aub_Pr_mean_cov",
                 "H2Aub_Pr_EZH2i_mean_cov")

    ub_colors <- colorRampPalette(c("white", "#8140c1"))(n = 1000)

    ub_subgroup <- values[indices, ub_cols]

    rownames(ub_subgroup) <- ub_subgroup$name
    ub_subgroup$name <- NULL
    colnames(ub_subgroup) <- c("H2Aub_Ni" ,
                               "H2Aub_Ni_EZH2i",
                               "H2Aub_Pr",
                               "H2Aub_Pr_EZH2i")

    ub_labels <- ub_subgroup[rev(genes_ordered), ]


    ub <- heatmaply(ub_subgroup[rev(genes_ordered), ], Rowv = F, Colv = F,
                    colors = ub_colors,
                    fontsize_row = row_font, fontsize_col = col_font,
                    colorbar_xanchor='left',colorbar_yanchor='top', colorbar_xpos=1.1, colorbar_ypos=0, colorbar_thickness = 30, plot_method = "plotly",
                    key.title = "H2Aub", limits = ub_limits)



    fig <- subplot(ub, k4, k27, h, widths = c(0.17, 0.17, 0.17, 0.49), shareY = T)
    fig
  }
