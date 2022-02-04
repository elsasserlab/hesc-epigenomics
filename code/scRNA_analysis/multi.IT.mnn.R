#' ---
#' title: "mnn integration"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---



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
# suppressMessages(library(scran))
# suppressMessages(library(batchelor))
suppressMessages(library(SeuratWrappers))

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

# working directory
DIR <- "/home/chenzh/My_project/Nerges_EZH2i_ScRNA"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
knitr::opts_chunk$set(echo=FALSE)


TD="Nov5_2021"

nGene=2500; pc=30

# Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
rename <- dplyr::rename

options(digits = 4)
options(future.globals.maxSize= 3001289600)

counts.filter <- readRDS(paste0("tmp_data/",TD,"/combineRefH.counts.filter.rds"))
meta.filter <- readRDS(paste0("tmp_data/",TD,"/combineRefH.meta.filter.rds"))%>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/combineRefH.lognormExp.mBN.rds"))



#' mnn integration
temp.M <- meta.filter 
temp.sel.expG <-rownames(lognormExp.mBN )

data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
#' release memory
rm(lognormExp.mBN,data.merge,counts.filter)

data.spt <- data.spt[c("JPF2019","D3post","SPH2016","CS7","Nerges_EZH2i_Naive","Nerges_EZH2i_Primed")]
set.seed(123)

# pdf("temp1.pdf")
# for (nGene in c(2000,2500,3000,3500,4000,4500,5000)) {
#   data.ob <- RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F)
#    print(nGene)
#   AP(data.ob)
#   APPJ(data.ob,"IBD2")
# }
# dev.off()
data.ob <- RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc)



data.temp <- data.ob %>% FindClusters(reso=0.4,verbose=F)

data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:mt.perc)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))%>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")

#' using old sph2016 annotation
cell_lineage_data <- read.delim("~/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.short.txt",stringsAsFactors = F,row.names = 1)%>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% setNames(c("cell","Embryo","Days","Lineage1","Lineage2")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="primitive endoderm", "PrE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="epiblast", "Epiblast")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="trophectoderm", "TE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="not applicable", "Prelineage")) %>% mutate(EML=Lineage1)%>% select(cell,EML) %>% filter(cell %in% data.ob.umap$cell)

data.ob.umap <- data.ob.umap %>% rows_update(cell_lineage_data,by="cell")

#'  format EML annotation
data.ob.umap <- data.ob.umap %>% mutate(EML=ifelse(EML=="EM_NA","Prelineage",EML)) %>% mutate(EML=ifelse(EML=="EPI","Epiblast",EML))%>% mutate(EML=ifelse(EML %in% c("PE","PrE"),"Endoderm",EML))  %>% mutate(EML=ifelse(pj=="D3post" & EML=="ICM","3D_ICM",EML))

#' rename EML annotation
data.ob.umap <- data.ob.umap %>% mutate(rename_EML=EML) %>% mutate(rename_EML=ifelse(rename_EML%in% c("AdvMes","AxMes","EmMes","NasMes"),"Mesoderm",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("CTB","EVT","STB","TE","EarlyTE"),"TE",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("MeLC2","MeLC1"),"MeLC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("Tsw-AMLC","AMLC"),"AMLC",rename_EML))

#' create cluster EML annotation
data.ob.umap <- data.ob.umap %>% filter(pj %in% c("Nerges_EZH2i_Naive")) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C0","C2","C4","C12"),"ELC","Undef")) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C9"),"HLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C7","C11"),"TLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C1"),"MeLC",cluster_EML)) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C10"),"MeLC",cluster_EML))  %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C3","C5"),"AMLC",cluster_EML)) %>% bind_rows(data.ob.umap %>% filter(!pj %in% c("Nerges_EZH2i_Naive")) %>% mutate(cluster_EML=rename_EML))


saveRDS(data.ob ,file=paste0("tmp_data/",TD,"/multi.IT.mnn.rds"))
saveRDS(data.ob.umap ,file=paste0("tmp_data/",TD,"/multi.IT.umap.cord.rds"))



#' some check
check=FALSE
if (check) {
  data.ob.umap  %>% filter(pj=="Nerges_EZH2i_Naive") %>% group_by(pj,EML,cluster_EML) %>% summarise(nCell=n()) %>%  ungroup() %>% select(EML,cluster_EML,nCell)
  temp <- data.ob.umap  %>% filter(pj=="Nerges_EZH2i_Naive") %>% mutate(cla=EML) %>% group_by(pj,devTime,cla) %>% summarise(nCell=n()) %>% group_by(pj,devTime) %>% mutate(prop=nCell/sum(nCell)) %>% mutate(devTime=gsub("EZH2i_Naive_","",devTime)) %>% mutate(devTime=ifelse(devTime=="WT","D0",devTime)) %>% mutate(devTime=factor(devTime,c("D0","D2","D4","D7"),ordered = T)) %>% ungroup() %>% arrange(devTime)
  print(
    temp %>% ggplot+geom_bar(mapping=aes(x=devTime,fill=cla,y=prop*100),stat="identity",position="dodge")+ theme_classic() +xlab("")+ylab("% cells")+ggtitle("Nerges_EZH2i_Naive")+FunTitle()
  )
  temp %>% select(devTime,cla,nCell)
  
  temp <- data.ob.umap  %>% filter(pj=="Nerges_EZH2i_Naive") %>% mutate(cla=cluster_EML) %>% group_by(pj,devTime,cla) %>% summarise(nCell=n()) %>% group_by(pj,devTime) %>% mutate(prop=nCell/sum(nCell)) %>% mutate(devTime=gsub("EZH2i_Naive_","",devTime)) %>% mutate(devTime=ifelse(devTime=="WT","D0",devTime)) %>% mutate(devTime=factor(devTime,c("D0","D2","D4","D7"),ordered = T)) %>% ungroup() %>% arrange(devTime)
  print(
    temp %>% ggplot+geom_bar(mapping=aes(x=devTime,fill=cla,y=prop*100),stat="identity",position="dodge")+ theme_classic() +xlab("")+ylab("% cells")+ggtitle("Nerges_EZH2i_Naive")+FunTitle()
  )
  temp %>% select(devTime,cla,nCell)
  
  DefaultAssay(data.temp) <- "RNA"
  
  VlnPlot(subset(data.temp,cells=(meta.filter %>% filter(pj=="Nerges_EZH2i_Naive")) %>% pull(cell)),c("ISL1","GABRP","IGFBP2"),ncol=1)
  VlnPlot(subset(data.temp,cells=(meta.filter %>% filter(pj=="JPF2019")) %>% pull(cell)),c("ISL1","GABRP","IGFBP2"),ncol=1)
  
  
  AP(data.ob)
  APPJ(data.ob,"SPH2016")
  APPJ(data.ob,"D3post")
  APPJ(data.ob,"CS7")
  APPJ(data.ob,"Nerges_EZH2i_Naive")
  APPJ(data.ob,"Nerges_EZH2i_Primed")
}



