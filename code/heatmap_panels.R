
library(wigglescout)
library(rtracklayer)
library(cowplot)
library(ggplot2)
#
# params <- list(datadir = "~/work/publications/hesc-epigenomics/data/")
# bed1 <- "~/work/publications/hesc-epigenomics/output/Kumar_2020/tss_2kb_up_k27m3_primed.bed"
# bed2 <- "~/work/publications/hesc-epigenomics/output/Kumar_2020/tss_2kb_down_k27m3_primed.bed"
#
# genes_signif_up <- import("~/work/publications/hesc-epigenomics/output/Kumar_2020/tss_2kb_up_k27m3_primed.bed")
# genes_signif_down <- import("~/work/publications/hesc-epigenomics/output/Kumar_2020/tss_2kb_down_k27m3_primed.bed")
#
#
# bwdir <- file.path(params$datadir, "bw/Kumar_2020")
# bwfiles <- list(k4 = list.files(bwdir, pattern = "H3K4m3.*pooled.hg38.scaled.*", full.names = T),
#                 k27 = list.files(bwdir, pattern = "H3K27m3.*pooled.hg38.scaled.*", full.names = T),
#                 ub = list.files(bwdir, pattern = "H2Aub.*pooled.hg38.scaled.*", full.names = T),
#                 input = list.files(bwdir, pattern = "IN.*pooled.hg38.*", full.names = T))
#
#
# bwlist <- bwfiles$k4

plot_bw_profile_row <- function(bwfiles, loci, bwnames, locinames, global_scale = FALSE, zmin = NULL, zmax = NULL, ...) {
  plot_list <- list()

  if (global_scale == TRUE) {
    limits <- calculate_profile_limits(bwfiles, loci, ...)
  }

  for (i in 1:length(bwfiles)) {
    plot_list[[i]] <- plot_bw_profile(bwfiles[[i]], loci, labels = locinames, verbose = FALSE, ...) + labs(title = "", x = "", y = "") + theme(plot.margin = unit(c(1,8,1,1), "mm"))
    if (global_scale == TRUE) {
      plot_list[[i]] <- plot_list[[i]] + ylim(limits[1], limits[2])
    }
    if (i < length(bwfiles)) {
      plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
    }
  }

  plot_list[[1]] <- plot_list[[1]] + labs(y = "RPGC")

  plot_list
}

plot_bw_heatmap_row <- function(bwfiles, loci, bwnames, lociname, global_scale = FALSE, zmin = NULL, zmax = NULL, ...) {
  ref_bw <- bwfiles[1] # or whichever is the reference

  h1 <- bw_heatmap(ref_bw, loci=loci, ...)
  ord <- order(rowMeans(h1[[1]]),decreasing = F)

  plot_list <- list()

  plot_list[[1]] <-
    plot_bw_heatmap(
      bwfiles[1],
      loci = loci,
      order_by = ord,
      verbose = FALSE,
      zmin = zmin,
      zmax = zmax,
      ...
    ) + labs(x = "", y = lociname, title = "") + theme(plot.margin = unit(c(1,8,1,1), "mm"))

  if (! is.null(bwnames)) {
    plot_list[[1]] <- plot_list[[1]] + labs(x = bwnames[1])
  }

  if (global_scale == TRUE) {
    plot_list[[1]] <- plot_list[[1]] + theme(legend.position = "none")
  }

  for (i in 2:length(bwfiles)) {
    plot_list[[i]] <-
      plot_bw_heatmap(
        bwfiles[i],
        loci = loci,
        order_by = ord,
        verbose = FALSE,
        zmin = zmin,
        zmax = zmax,
        ...
      ) + labs(x = "", y = "", title = "") + theme(plot.margin = unit(c(1,8,1,1), "mm"))
    if (!is.null(bwnames)) {
      plot_list[[i]] <- plot_list[[i]] + labs(x = bwnames[i])
    }
    if (global_scale == TRUE) {
      plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
    }
  }
  plot_list
}

calculate_matrix_limits <- function(bwfiles, lociset, ...) {
  mats <- list()
  for (i in 1:length(lociset)) {
    mats <- c(mats, bw_heatmap(bwfiles, loci = lociset[[i]], ...))
  }
  global_limits(mats)
}

calculate_profile_limits <- function(bwfiles, lociset, ...) {
  values <- c()
  for (i in 1:length(lociset)) {
    new_values <- bw_profile(bwfiles, loci = lociset[[i]], ...)
    values <- c(values, new_values$mean)
  }
  c(min(values, na.rm = TRUE), max(values, na.rm = TRUE))
}


global_limits <- function(matlist) {
  minima <- lapply(matlist, single_mat_lower_limits)
  maxima <- lapply(matlist, single_mat_upper_limits)
  c(min(unlist(minima)), max(unlist(maxima)))
}

single_mat_lower_limits <- function(m) {
  quantile(unlist(m), c(0.01), na.rm = TRUE)
}
single_mat_upper_limits <- function(m) {
  quantile(unlist(m), c(0.99), na.rm = TRUE)
}


# first bwfile is the reference
plot_bw_heatmap_panel <- function(bwfiles, lociset, bwnames, locinames, global_scale = FALSE, proportional = FALSE, zmin = NULL, zmax = NULL, ...) {

  profiles <- plot_bw_profile_row(bwfiles, loci = lociset, locinames = locinames, global_scale = global_scale, ...)

  if (!is.null(zmin) && !is.null(zmax)) {
    global_limits <- c(zmin, zmax)
  }
  else {
    if (global_scale == TRUE) {
      global_limits <- calculate_matrix_limits(bwfiles, lociset, ...)
      zmin <- global_limits[1]
      zmax <- global_limits[2]
    }
  }


  rows <- list()
  for (i in 1:length(lociset)) {
    row_names <- NULL
    if (i == length(lociset)) {
        row_names <- bwnames
    }
    new_row <- plot_bw_heatmap_row(bwfiles, lociset[[i]], row_names, locinames[i], global_scale = global_scale, zmin = zmin, zmax = zmax, ...)
    rows <- c(rows, new_row)
  }

  heatmap_rel_heights <- rep(0.4, length(lociset))

  if (proportional == TRUE) {
    min_possible_height <- 0.1
    heights <- sapply(lociset, length)
    heatmap_rel_heights <- heights / max(heights)
    heatmap_rel_heights <- heatmap_rel_heights * 0.4

    heatmap_rel_heights[heatmap_rel_heights < min_possible_height] <- min_possible_height

  }

  length_data <- data.frame(loci=locinames, size=heights)
  lengths_plot <- ggplot(length_data, aes(x=loci, y=size)) + geom_bar(stat = "identity") + theme_default(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", title = "Group size")

  panel <- plot_grid(plotlist = c(profiles, rows), nrow = length(lociset)+1, rel_heights = c(0.2, heatmap_rel_heights))

  if (global_scale == TRUE) {
    legend <- get_legend(rows[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12), legend.position = "right"))
    side_panel <- plot_grid(lengths_plot, legend, ncol = 1, rel_heights = c(0.2, sum(heatmap_rel_heights)))
  } else {
    side_panel <- lengths_plot
  }
  plot_grid(panel, side_panel, rel_widths = c(1, .25))

}

plot_lengths <- function(lociset, locinames) {
  length_data <- data.frame(loci=locinames, size=sapply(lociset, length))
  length_data$loci <- factor(length_data$loci, levels = locinames)
  lengths_plot <- ggplot(length_data, aes(x=loci, y=size)) +
    geom_bar(stat = "identity", color = "black") +
    theme_default(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x="", title = "Group size")
  lengths_plot
}
