merge_granges_plus_name <- function(granges_c1, granges_c2) {
  # Bind first, get numbers after
  granges_c1 <- sortSeqlevels(granges_c1)
  granges_c1 <- sort(granges_c1)

  granges_c2 <- sortSeqlevels(granges_c2)
  granges_c2 <- sort(granges_c2)

  fields <- names(mcols(granges_c1))

  if ("name" %in% fields) {
    granges_c1 <- granges_c1[, fields[fields != "name"]]
  }

  fields <- names(mcols(granges_c2))
  if ("name" %in% fields) {
    granges_c2 <- granges_c2[, fields[fields != "name"]]
  }

  cts_df <- cbind(data.frame(granges_c1), mcols(granges_c2))

  cts_df


}

get_names_values <- function(granges) {
  names_values <- NULL
  if ("name" %in% names(mcols(granges))) {
    names_values <- mcols(granges)[["name"]]
  }
  names_values
}

#' Get DESeq2 results for a pair of GRanges lists.
#'
#' Loci in condition 1 and 2 must match. If a name file exists, it will be kept.
#' This also needs to match. Names from condition 1 are kept.
#'
#' Non-complete cases are also dropped before calculation.
#'
#' @param granges_c1 GRanges values for condition 1
#' @param granges_c2 GRanges values for condition 2
#' @param label_c1 Name of condition 1
#' @param label_c2 Name of condition 2
#' @param estimate_size_factors Whether to estimate size factors or not. If
#'   TRUE, the behavior is the default on DESeq2 documentation. If FALSE, this
#'   step is skipped and all the factors are given weight = 1. This is done to
#'   provide support for quantitative Minute-ChIP data.
#' @param length_factor Scaling factor for read counts.
#' @param as_granges Whether to return a GRanges object.
#'
#' @return If as_granges == TRUE, a GRanges object. Otherwise, a DESeq results.
bw_granges_diff_analysis <- function(granges_c1,
                                     granges_c2,
                                     label_c1,
                                     label_c2,
                                     estimate_size_factors,
                                     length_factor = 100,
                                     as_granges = FALSE) {

  used_values <- function(gr, skip_fields = c("name", "gene_id")) {
    names_list <- names(mcols(gr))
    length(setdiff(names_list, skip_fields))
  }

  cts_df <- merge_granges_plus_name(granges_c1, granges_c2)
  names_values <- get_names_values(granges_c1)
  if (! is.null(names_values)) {
    rownames(cts_df) <- names_values
  }
  # Needs to drop non-complete cases and match rows
  complete <- complete.cases(cts_df)
  cts_df <- cts_df[complete, ]

  values_df <- cts_df[, 6:ncol(cts_df)] %>% dplyr::select(where(is.numeric))
  cts <- get_nreads_columns(values_df, cts_df$width, factor = length_factor)

  condition_labels <- c(rep(label_c1, used_values(granges_c1)),
                        rep(label_c2, used_values(granges_c2)))

  coldata <- data.frame(colnames(cts), condition = as.factor(condition_labels))
  coldata$condition <- relevel(coldata$condition, ref = label_c1)

  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition,
                                rowRanges = granges_c1[complete, ])

  if (estimate_size_factors == TRUE) {
    dds <- estimateSizeFactors(dds)
  }
  else {
    sizeFactors(dds) <- c(rep(1, ncol(cts)))
  }

  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)

  if (as_granges) {
    result <- results(dds, format = "GRanges")
    if (!is.null(names_values)) {
      result$name <- names_values[complete]
    }

  }
  else {
    result <- dds
  }

  result
}

get_nreads_columns <- function(df, lengths, fraglen = 150, factor = 1) {
  # Convert mean coverages to round integer read numbers
  cts <- as.matrix(df)

  # Estimate of # fragments from coverage
  cts <- round(cts*factor*(lengths / fraglen))
  cts
}




rsem_deseq_analysis <- function(counts_file, c1_columns, c2_columns, c1_name, c2_name, reference, alpha = 0.05, shrink = F) {
  if (! reference %in% c(c1_name, c2_name)) {
    stop(paste(reference, "must be", c1_name, "or", c2_name))
  }

  counts <- read.table(counts_file, sep = "\t", header = T)

  columns <- c(c1_columns, c2_columns)

  samples <- data.frame(row.names = columns, condition = factor(c(rep(c1_name, length(c1_columns)), rep(c2_name, length(c2_columns)))))
  samples$condition <- relevel(samples$condition, ref = reference)

  counts_only <- round(counts[, columns])
  rownames(counts_only) <- counts$gene_id

  dds <- DESeqDataSetFromMatrix(countData = counts_only,
                                colData = samples,
                                design = ~ condition)

  coef_name <- paste("condition", setdiff(c(c1_name, c2_name), reference), "vs", reference, sep = "_")

  deseq_analysis(dds, coef_name, alpha = alpha, shrink = shrink)
}

deseq_condition_table <- function(c1_columns, c2_columns, c1_name, c2_name, reference) {
  columns <- c(c1_columns, c2_columns)

  samples <- data.frame(row.names = columns, condition = factor(c(rep(c1_name, length(c1_columns)), rep(c2_name, length(c2_columns)))))
  samples$condition <- relevel(samples$condition, ref = reference)
  samples
}

deseq_analysis <- function(dds, coef_name, alpha = 0.05, shrink = F) {
  dds <- DESeq(dds)
  res <- NULL
  if (shrink == FALSE) {
    res <- results(dds, alpha=alpha)
  } else {
    res <- lfcShrink(dds, coef=coef_name, type="apeglm")
  }

  res
}


