#' ---
#' title: "merge other reference dataset and scran norm"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---




# ### Loading R library

#R3.6
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(scran))
suppressMessages(library(batchelor))

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

# working directory
DIR <- "/home/chenzh/My_project/Ploycomb"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
knitr::opts_chunk$set(echo=FALSE)


TD="Nov5_2021"

# Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
rename <- dplyr::rename

options(digits = 4)
options(future.globals.maxSize= 3001289600)


#' loading human reference dataset
#' please check https://github.com/zhaocheng3326/CheckBlastoids_scripts for other reference 
meta.ref.filter <- readRDS("/home/chenzh/My_project/JP_project/tmp_data/Nov14_2021/meta.filter.rds") %>% filter(pj %in% c("SPH2016","D3post","CS7","JPF2019"))
counts.ref.filter <- readRDS("/home/chenzh/My_project/JP_project/tmp_data/Nov14_2021/counts.filter.rds")[,meta.ref.filter$cell]


#' loading EZH2i dataset
EZH2i.meta <- readRDS(paste0("tmp_data/",TD,"/meta.filter.EZH2i.rds"))
EZH2i.counts <-  readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))

#' check
colnames(EZH2i.meta) %>% setdiff(colnames(meta.ref.filter ))
intersect(EZH2i.meta$cell,meta.ref.filter$cell)
(rownames(EZH2i.counts)==rownames(counts.ref.filter)) %>% table()

#' merge
counts.filter <- counts.ref.filter %>% bind_cols(EZH2i.counts ) %>% as.data.frame()
meta.filter <- meta.ref.filter %>% bind_rows(EZH2i.meta)

#' normalization
#' get the expressed genes ( expressed in at least 1 datatsets)
expG.set <- list()
for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
  temp.cell <- meta.filter %>% filter(pj==b) %>% pull(cell)
  expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
}
sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()

sce.ob <- list()
for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
  print(b)
  temp.M <- meta.filter %>% filter(pj==b) 
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
  sce.ob[[b]] <- temp.sce
}

mBN.sce.ob <- multiBatchNorm(sce.ob$Nerges_EZH2i_Naive,sce.ob$Nerges_EZH2i_Primed,sce.ob$JPF2019,sce.ob$SPH2016,sce.ob$D3post,sce.ob$CS7)
lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)


saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/combineRefH.counts.filter.rds"))
saveRDS(meta.filter,file=paste0("tmp_data/",TD,"/combineRefH.meta.filter.rds"))
saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/combineRefH.lognormExp.mBN.rds"))