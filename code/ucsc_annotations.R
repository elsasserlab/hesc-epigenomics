#' Get GRanges with gene symbols for UCSC hg38
#'
#' @return A GRanges where name field is Gene symbol
genes_hg38 <- function() {
  # This matches gene IDs to symbols
  gene_symbol_db <- org.Hs.egSYMBOL

  # Get the gene names that are mapped to an entrez gene identifier
  mapped_genes <- mappedkeys(gene_symbol_db)
  gene_info <- as.list(gene_symbol_db[mapped_genes])

  # Apparently there are some pathway things that put NAs in the list
  gene_info <- gene_info[!is.na(gene_info)]
  gene_id_to_name <- unlist(gene_info)

  genes <- suppressMessages(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  genes$name <- unname(gene_id_to_name[genes$gene_id])

  sort(genes[!is.na(genes$name),])
}

#' Get UCSC hg38 genes by gene symbol
#'
#' @param gene_names List of gene symbols
#'
#' @return A GRanges object with gene_id and name metadata columns
get_genes_by_name <- function(gene_names) {
  all_genes <- genes_hg38()
  all_genes[all_genes$name %in% gene_names, ]
}

#' Get UCSC hg38 genes that overlap with a set of loci
#'
#' @param loci GRanges object
#' @param minoverlap Minimum overlap (0 by default).
#'
#' @return A GRanges object with gene_id and name metadata columns
get_genes_by_overlap <- function(loci, minoverlap = 0L) {
  all_genes <- genes_hg38()
  subsetByOverlaps(all_genes, loci, minoverlap = minoverlap)
}

