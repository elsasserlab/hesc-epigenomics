#' ---
#' title: "QC"
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
DIR <- "/home/chenzh/My_project/Nerges_EZH2i_ScRNA"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
knitr::opts_chunk$set(echo=FALSE)


# Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
rename <- dplyr::rename

options(digits = 4)
options(future.globals.maxSize= 3001289600)


qc.nGene.min <- 1000
qc.nGene.max <- 6000
qc.mt.perc <- 0.15


TD="Nov5_2021"


if (file.exists(paste0("tmp_data/",TD,"/meta.filter.rds"))) {
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
}else{
  #' loading data
  
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  load(paste0("tmp_data/",TD,"/gene.meta.Rdata"),verbose=T)
  
  
  
  meta.filter <- meta.all %>% filter( nGene < qc.nGene.max & mt.perc < qc.mt.perc & nGene > qc.nGene.min) %>% mutate_all(as.vector) %>% mutate(pj=ifelse(pj=="Nerges_EZH2i_NaiveD7","Nerges_EZH2i_Naive",pj))
  
  counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene),meta.filter$cell]
  
  
  #' get the expressed genes ( expressed in at least 1 datatsets)
  # expG.set <- list()
  # for (b in unique(meta.filter$devTime  %>% unique() %>% as.vector())) { 
  #   temp.cell <- meta.filter %>% filter(devTime==b) %>% pull(cell)
  #   expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  # }
  # sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  # 
  # sce.ob <- list()
  # for (b in unique(meta.filter$devTime  %>% unique() %>% as.vector())) { 
  #   print(b)
  #   temp.M <- meta.filter %>% filter(devTime==b) 
  #   temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
  #   sce.ob[[b]] <- temp.sce
  # }
  # 
  # mBN.sce.ob <- multiBatchNorm(sce.ob$EZH2i_Naive_D2,sce.ob$EZH2i_Naive_D4,sce.ob$EZH2i_Naive_WT,sce.ob$EZH2i_Primed_D7,sce.ob$EZH2i_Primed_WT,sce.ob$EZH2i_Naive_D7)
  # lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  # 
  
  saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/counts.filter.rds"))
  saveRDS(meta.filter,file=paste0("tmp_data/",TD,"/meta.filter.rds"))
  #saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
 
}


if (file.exists(paste0("tmp_data/",TD,"/lognormExp.mBN.NPonly.rds"))) {
  print("done")
}else{
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
 
  expG.set <- list()
  for (b in c("Nerges_EZH2i_Naive","Nerges_EZH2i_Primed")) { 
    temp.cell <- meta.filter %>% filter(pj==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  
  
  sce.ob <- list()
  for (b in c("Nerges_EZH2i_Naive","Nerges_EZH2i_Primed")) { 
    print(b)
    temp.M <- meta.filter %>% filter(pj==b) 
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
    sce.ob[[b]] <- temp.sce
  }
  
  mBN.sce.ob <- multiBatchNorm(sce.ob$Nerges_EZH2i_Naive,sce.ob$Nerges_EZH2i_Primed)
  lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/lognormExp.mBN.NPonly.rds"))
}


#' #### check the general distribution
meta.all %>% ggplot+geom_histogram(mapping=aes(x=nGene,fill=devTime),bins=100)+geom_vline(xintercept=qc.nGene.min)+geom_vline(xintercept=qc.nGene.max)+facet_wrap(.~devTime)+ylab("Number of cells")
meta.all %>% ggplot+geom_boxplot(mapping=aes(x=devTime,y=nGene,fill=devTime))+ theme(axis.text.x=element_text(angle = 90))
meta.all %>% ggplot+geom_histogram(mapping=aes(x=mt.perc,fill=devTime),bins=100)+geom_vline(xintercept=qc.mt.perc)+facet_wrap(.~devTime)+ylab("Number of cells")
meta.all %>% ggplot+geom_point(mapping=aes(x=nGene,y=mt.perc,col=devTime))+facet_wrap(.~devTime)


#' #### number of cells
#+ fig.width=9,fig.height=9
meta.all %>% group_by(devTime) %>% summarise(BeforeQCnCell=n()) %>% inner_join(meta.filter  %>% group_by(devTime) %>% summarise(AfterQCnCell=n()),by="devTime") %>% gather(QC,nCell,-devTime) %>% ggplot+geom_bar(mapping=aes(x=devTime,fill=QC,y=nCell),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90))

#' check the cell number during QC
table(meta.filter$devTime)
table(meta.all$devTime)