library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(biomaRt)

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

trx_ensembl_to_refseq <- function(gene_id) {
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  result <- getBM(attributes=c("refseq_mrna", "ensembl_transcript_id", "hgnc_symbol"), filters = "ensembl_transcript_id", values = gene_id, mart= ensembl)
  result
}

canonical_genes_hg38 <- function(file_ref_flat, file_canonical) {
  # These files are in the default format I got
  # refFlat.txt from UCSC annotation:
  # knownCanonical.txt from Ensembl
  ref_flat <- read.table(file_ref_flat, sep = "\t")
  canonical <- read.table(file_canonical, sep = "\t")

  # remove the version
  canonical_ids <- gsub("[[:punct:]][0-9]+", "", canonical[, 5])

  refseq_ids <- trx_ensembl_to_refseq(canonical_ids)

  ref_flat <- ref_flat[, c(1:6)]
  columns <- c("name", "id", "chr", "strand", "start", "end")
  colnames(ref_flat) <- columns


  canonical_flat <- ref_flat[ref_flat[, "id"] %in% refseq_ids$refseq_mrna, ]

  # I still seem to need to select one each, because in ensembl_to_refseq
  # there are sometimes 1:n correspondences.
  canonical_flat$length <- canonical_flat$end - canonical_flat$start
  canonical_flat <- canonical_flat %>% group_by(name) %>% slice_max(length, n = 1, with_ties = F)

  canonical_genes <- unique(canonical_flat$name)

  not_genes <- setdiff(unique(ref_flat[, "name"]), canonical_genes)

  non_canonical_genes <- ref_flat %>% filter(!name %in% canonical_genes)
  non_canonical_genes$length <- non_canonical_genes$end - non_canonical_genes$start
  non_canonical_genes <- non_canonical_genes %>% group_by(name) %>% slice_max(length, n = 1, with_ties = F)

  all <- rbind(canonical_flat[, columns], non_canonical_genes[, columns]) %>% arrange(name)
  all$id <- NULL
  makeGRangesFromDataFrame(all, keep.extra.columns = T)
}

genes_from_rnaseq <- function(refgene) {
  table <- read.table(refgene, sep = "\t")
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

