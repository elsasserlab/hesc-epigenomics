# R3.6
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

#' workding directory
argv=c("/home/chenzh/My_project/Nerges_EZH2i_ScRNA")
setwd(argv[1])
knitr::opts_knit$set(root.dir=argv[1])
knitr::opts_chunk$set(echo=FALSE)
setwd(argv[1])


# working directory

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(monocle))
#suppressMessages(library(monocle3))
suppressMessages(library(slingshot))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
#source("loca.quick.fun.R")

TD="Nov5_2021"
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.EZH2i.rds"))
counts.filter <-  readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/multi.IT.umap.cord.rds"))%>% mutate(cluster_EML=ifelse(cluster_EML=="ExE_MeLC" & pj=="Nerges_EZH2i_Naive","MeLC",cluster_EML))

Gene.discrp<- read.delim("~/Genome/Human/RefSeq/Homo_sapiens.GRCh38.gene.description",sep="\t",stringsAsFactors=F,header = F,row.names=1) %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% rename(anno=V2)
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.NPonly.rds"))

options(digits = 4)
options(future.globals.maxSize= 3001289600)

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

if (file.exists(paste0("tmp_data/",TD,"/naive.trajectory.Rdata"))) {
  load(paste0("tmp_data/",TD,"/naive.trajectory.Rdata"),verbose = T)
}else{
  
  temp.M <-  data.ob.umap %>% filter( pj=="Nerges_EZH2i_Naive")
  temp.sel.expG <- rownames(counts.filter )[rowSums(counts.filter[, temp.M$cell] >=1) >=5]
  
  data.temp <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 1000, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:25,verbose=F) %>% FindNeighbors( dims = 1:25,verbose = FALSE)
  
  cds <- newimport(data.temp) #### transfer the data to monocle object for monocle2
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  cds.diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~cluster_EML",cores=30) %>% arrange(qval)
  disp_table <- dispersionTable(cds)
  #cds.traj.gene <- disp_table %>% tbl_df() %>% filter(mean_expression >= 0.1 & dispersion_empirical > 1 * dispersion_fit) %>% pull(gene_id)
  cds.traj.gene <- row.names (subset(cds.diff_test_res %>% arrange(qval), qval < 10^(-20)))  
  length( cds.traj.gene )
  cds.out <- setOrderingFilter(cds, cds.traj.gene  )
  cds.out <- reduceDimension(cds.out, max_components = 2,reduction_method = 'DDRTree', verbose = F)
  cds.out <- orderCells(cds.out)
  
  #cds.out.psd.diff_test_res <- differentialGeneTest(cds.out,fullModelFormulaStr = "~sm.ns(Pseudotime)")
  cds.out.DM <- t(cds.out@reducedDimS) %>% as.data.frame()%>% tbl_df() %>% mutate(cell=colnames(cds.out)) %>% rename(traj_Dim1=V1,traj_Dim2=V2) %>% left_join(pData(cds.out) %>% tibble::rownames_to_column("cell") %>% select(cell,SID:State) %>% mutate(State=as.vector(State)) %>% tbl_df() ,by="cell") 
  
  psd.cutoff.topN <- ((cds.out.DM %>% filter(cluster_EML=="ELC") %>% filter(devTime=="EZH2i_Naive_WT") %>% nrow()) * 0.05) %>% round()
  psd.cutoff <- cds.out.DM %>% filter(cluster_EML=="ELC") %>% filter(devTime=="EZH2i_Naive_WT") %>% arrange(desc(Pseudotime)) %>% head(psd.cutoff.topN) %>% tail(1) %>% pull(Pseudotime)
  
  
  cds.out.DM <- cds.out.DM %>% mutate(cluster_EML_state=cluster_EML) %>% mutate(cluster_EML_state=ifelse(Pseudotime < psd.cutoff & cluster_EML=="ELC", "gELC",cluster_EML_state))%>% mutate(cluster_EML_state=ifelse(Pseudotime >= psd.cutoff & cluster_EML=="ELC" & State=="1" ,   "aELC" ,cluster_EML_state)) %>% mutate(cluster_EML_state=ifelse(Pseudotime >= psd.cutoff & cluster_EML=="ELC" & State=="2","TaELC" , cluster_EML_state)) %>% mutate(cluster_EML_state=ifelse(Pseudotime >= psd.cutoff & cluster_EML=="ELC" & State=="3","MaELC" , cluster_EML_state)) 
  pData(cds.out)$cluster_EML_state <- (cds.out.DM %>% tibble::column_to_rownames("cell"))[rownames(pData(cds.out)),"cluster_EML_state"]
  psd.branch <- cds.out.DM %>% filter(cluster_EML_state=="aELC")%>% pull(Pseudotime) %>% max()
  
  
  #' create psd-pseduo_time Seurat_object
  data.psd.ob <- data.temp
  data.psd.ob@reductions$umap@cell.embeddings[,"UMAP_1"] <- (cds.out.DM %>% tibble::column_to_rownames("cell"))[rownames(data.psd.ob@meta.data),"traj_Dim1"]
  data.psd.ob@reductions$umap@cell.embeddings[,"UMAP_2"] <- (cds.out.DM %>% tibble::column_to_rownames("cell"))[rownames(data.psd.ob@meta.data),"traj_Dim2"]
  data.psd.ob@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.psd.ob)])
  
  #' save object
  save(psd.cutoff,psd.branch,file=paste0("tmp_data/",TD,"/psd.cutoff.Rdata"))
  saveRDS(cds.out.DM, file=paste0("tmp_data/",TD,"/naive.trajectory.dm.rds"))
  saveRDS(data.psd.ob,file=paste0("tmp_data/",TD,"/psd.seurat.psdTime.rds"))
  save(cds, cds.diff_test_res,  cds.out, disp_table, cds.out.DM,file=paste0("tmp_data/",TD,"/naive.trajectory.Rdata"))
}
load(paste0("tmp_data/",TD,"/psd.cutoff.Rdata"),verbose=T)


cowplot::plot_grid(
  plot_cell_trajectory(cds.out, color_by = "State",show_branch_points=F)+NoAxes(),
  plot_cell_trajectory(cds.out, color_by = "cluster_EML")+NoAxes(),
  plot_cell_trajectory(cds.out, color_by = "Pseudotime")+NoAxes(),
  plot_cell_trajectory(cds.out, color_by = "cluster_EML_state")+NoAxes()
)

plot_cell_trajectory(cds.out, color_by = "cluster_EML_state")+NoAxes()+geom_vline(xintercept = psd.cutoff ,linetype="dashed")

#plot_genes_branched_pseudotime(cds.out[c("GATA2"),],branch_point = 1,color_by = "cluster_EML", ncol = 1)

#' check the Pseudotime density
print(
  cds.out.DM %>% mutate(od=factor(devTime,c("EZH2i_Naive_WT","EZH2i_Naive_D2","EZH2i_Naive_D4","EZH2i_Naive_D7"),ordered = T)) %>% arrange(od) %>% mutate()%>% ggplot(mapping=aes(x=Pseudotime,y = od,fill=devTime))+ggridges::geom_density_ridges()+theme_classic()+ggtitle("EZH2i_Naive cells")+FunTitle()+ylab("")+NoLegend()#+facet_wrap(~State)
)






