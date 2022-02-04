# R_test
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R_test"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

#' workding directory
argv=c("/home/chenzh/My_project/Nerges_EZH2i_ScRNA")
setwd(argv[1])
knitr::opts_knit$set(root.dir=argv[1])
knitr::opts_chunk$set(echo=FALSE)
setwd(argv[1])

TD="Nov5_2021"


# working directory

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(pheatmap))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
#source("loca.quick.fun.R")

options(digits = 4)
options(future.globals.maxSize= 3001289600)

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter


#' loading R object
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.EZH2i.rds"))
counts.filter <-  readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.NPonly.rds"))
sel.expG <- rownames(lognormExp.mBN)
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/multi.IT.umap.cord.rds"))
Gene.discrp<- read.delim("~/Genome/Human/RefSeq/Homo_sapiens.GRCh38.gene.description",sep="\t",stringsAsFactors=F,header = F,row.names=1) %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% rename(anno=V2)
cds.out.DM <- readRDS(paste0("tmp_data/",TD,"/naive.trajectory.dm.rds"))
load(paste0("tmp_data/",TD,"/psd.cutoff.Rdata"),verbose=T)


#' create Naive WT and Primed WT only object
temp.M <- meta.filter %>% filter(devTime %in% c("EZH2i_Naive_WT","EZH2i_Primed_WT"))
temp.sel.expG <-rownames(lognormExp.mBN )
data.np <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.np@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.np)])
Idents(data.np) <- as.factor(data.np@meta.data$devTime)

#' create Naive WT and treated object 
temp.M <- cds.out.DM  %>% filter(pj %in% c("Nerges_EZH2i_Naive")) %>% filter(cluster_EML!="Undef")
temp.sel.expG <-rownames(lognormExp.mBN )
data.naive <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.naive@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.naive)])
Idents(data.naive) <- as.factor(data.naive@meta.data$cluster_EML_state)

#' temp removed quick check
tg <- c("CDX1","HAND1","MIXL1","TBX3","TBXT","GATA3")


#cell.list[["D7_Naive_TLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML=="TLC") %>% pull(cell)
# cell.list[["D7_Naive_MeLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML=="MeLC") %>% pull(cell)
#cell.list[["D4_Naive_ELC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML=="ELC") %>% pull(cell)
# cell.list[["D4_Naive_TLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML=="TLC") %>% pull(cell)
# cell.list[["D4_Naive_MeLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML=="MeLC") %>% pull(cell)
# cell.list[["D7_Naive"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" ) %>% pull(cell)
# cell.list[["D7_Primed"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Primed_D7" ) %>% pull(cell)
# cell.list[["WT_Naive"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_WT" ) %>% pull(cell)
# cell.list[["WT_Primed"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Primed_WT" ) %>% pull(cell)
#cell.list[["D7_TaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="TaELC") %>% pull(cell)
#cell.list[["D7_MaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="MaELC") %>% pull(cell)
#cell.list[["D4_TaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML_state=="TaELC") %>% pull(cell)
#cell.list[["D4_MaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML_state=="MaELC") %>% pull(cell)

#cell.list[["D7_gELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="gELC") %>% pull(cell)
#cell.list[["D2_gELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D2" & cluster_EML_state=="gELC") %>% pull(cell)
#cell.list[["WT_gELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_WT" & cluster_EML_state=="gELC") %>% pull(cell)
#cell.list[["D7_aELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="aELC") %>% pull(cell)
#cell.list[["D4_aELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML_state=="aELC") %>% pull(cell)

#cell.list[["aELC"]] <- cds.out.DM %>% filter( cluster_EML_state=="aELC") %>% pull(cell)
#cell.list[["gELC"]] <- cds.out.DM %>% filter( cluster_EML_state=="gELC") %>% pull(cell)

#cell.list[["TaELC"]] <- cds.out.DM %>% filter( cluster_EML_state=="TaELC") %>% pull(cell)
#cell.list[["MaELC"]] <- cds.out.DM %>% filter( cluster_EML_state=="MaELC") %>% pull(cell)
# cell.list[["D7_Naive_MeLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML=="MeLC") %>% pull(cell)
# cell.list[["D4_Naive_TLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML=="TLC") %>% pull(cell)
# cell.list[["D4_Naive_MeLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML=="MeLC") %>% pull(cell)
# cell.list[["D7_Naive"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" ) %>% pull(cell)
# cell.list[["D7_Primed"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Primed_D7" ) %>% pull(cell)
# cell.list[["WT_Naive"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_WT" ) %>% pull(cell)
# cell.list[["WT_Primed"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Primed_WT" ) %>% pull(cell)

#paste("D7_gELC","vs","WT_gELC",sep="_"),
#paste("D4_gELC","vs","WT_gELC",sep="_"),
#paste("aELC","vs","gELC",sep="_"),
#paste("MaELC","vs","aELC",sep="_"),
#paste("TaELC","vs","aELC",sep="_"),
#paste("TaELC","vs","MaELC",sep="_"),
#paste("Naive_TLC","vs","Naive_MeLC",sep="_"),
#paste("D7_aELC","vs","D7_gELC",sep="_"),
#paste("D4_aELC","vs","D4_gELC",sep="_")
#paste("D7_S2T_ELC","vs","D7_S3M_ELC",sep="_"),
#paste("D4_S2T_ELC","vs","D4_S3M_ELC",sep="_"),
#paste("D7_Naive_TLC","vs","D7_Naive_ELC",sep="_"),
#paste("D7_Naive_MeLC","vs","D7_Naive_ELC",sep="_"),
#paste("D4_Naive_TLC","vs","D4_Naive_ELC",sep="_"),
#paste("D4_Naive_MeLC","vs","D4_Naive_ELC",sep="_"),
#paste("D7_Naive_ELC","vs","WT_Naive",sep="_"),
#paste("D7_Naive_TLC","vs","WT_Naive",sep="_"),
#paste("D7_Naive_MeLC","vs","WT_Naive",sep="_"),
#paste("D7_Primed","vs","WT_Primed",sep="_"),
#paste("D7_Naive","vs","WT_Naive",sep="_")
#paste("D7_MaELC","vs","D7_aELC",sep="_"),
#paste("D4_TaELC","vs","D4_aELC",sep="_"),

cell.list <- list()
cell.list[["WT_gELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_WT" & cluster_EML_state=="gELC") %>% pull(cell)
cell.list[["D4_gELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML_state=="gELC") %>% pull(cell)
cell.list[["D4_aELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML_state=="aELC") %>% pull(cell)
cell.list[["D4_TaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML_state=="TaELC") %>% pull(cell)
cell.list[["D4_MaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML_state=="MaELC") %>% pull(cell)

cell.list[["D7_gELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="gELC") %>% pull(cell)
cell.list[["D7_aELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="aELC") %>% pull(cell)
cell.list[["D7_TaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="TaELC") %>% pull(cell)
cell.list[["D7_MaELC"]] <- cds.out.DM %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML_state=="MaELC") %>% pull(cell)

cell.list[["Naive_aELC"]] <- cds.out.DM %>% filter( cluster_EML_state=="aELC") %>% pull(cell)
cell.list[["Naive_gELC"]] <- cds.out.DM %>% filter( cluster_EML_state=="gELC") %>% pull(cell)

cell.list[["D7_TLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML=="TLC") %>% pull(cell)
cell.list[["D7_MeLC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML=="MeLC") %>% pull(cell)
cell.list[["Naive_TLC"]] <- data.ob.umap %>% filter( cluster_EML=="TLC" & pj=="Nerges_EZH2i_Naive") %>% pull(cell)
cell.list[["Naive_MeLC"]] <- data.ob.umap %>% filter( cluster_EML=="MeLC" & pj=="Nerges_EZH2i_Naive") %>% pull(cell)

cell.list[["WT_ELC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_WT" & cluster_EML=="ELC") %>% pull(cell)
cell.list[["D2_ELC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D2" & cluster_EML=="ELC") %>% pull(cell)
cell.list[["D4_ELC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D4" & cluster_EML=="ELC") %>% pull(cell)
cell.list[["D7_ELC"]] <- data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7" & cluster_EML=="ELC") %>% pull(cell)

cell.list[["D47_gELC"]] <- cds.out.DM %>% filter(devTime %in% c("EZH2i_Naive_D4","EZH2i_Naive_D7") & cluster_EML_state=="gELC") %>% pull(cell)
cell.list[["D47_aELC"]] <- cds.out.DM %>% filter(devTime  %in% c("EZH2i_Naive_D4","EZH2i_Naive_D7") & cluster_EML_state=="aELC") %>% pull(cell)



lapply(cell.list,length)



DEG.stat=matrix(nrow=2,ncol=14)
rownames(DEG.stat) <- c("up_regulated","down_regulated")
colnames(DEG.stat) <- c(
  paste("D7_gELC","vs","WT_gELC",sep="_"),
  paste("D4_gELC","vs","WT_gELC",sep="_"),
  paste("D7_aELC","vs","D7_gELC",sep="_"),
  paste("D7_TaELC","vs","D7_MaELC",sep="_"),
  paste("D4_aELC","vs","D4_gELC",sep="_"),
  paste("D4_TaELC","vs","D4_MaELC",sep="_"),
  
  paste("D7_aELC","vs","WT_gELC",sep="_"),
  paste("Naive_aELC","vs","Naive_gELC",sep="_"),
  paste("D47_aELC","vs","D47_gELC",sep="_"),
  
  paste("D7_TLC","vs","D7_MeLC",sep="_"),
  paste("Naive_TLC","vs","Naive_MeLC",sep="_"),
  
  paste("D2_ELC","vs","WT_ELC",sep="_"),
  paste("D4_ELC","vs","WT_ELC",sep="_"),
  paste("D7_ELC","vs","WT_ELC",sep="_")
)



DEG.results <- list()
savefile <- paste0("tmp_data/",TD,"/DEG_branch_MAST.Rdata")
if (file.exists(savefile)){
  load(savefile,verbose = T)
}else{
  temp.M <- data.ob.umap %>% filter(pj %in% c("Nerges_EZH2i_Naive","Nerges_EZH2i_Primed"))
  temp.sel.expG <-rownames(lognormExp.mBN )
  
  data.ob.deg <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.ob.deg@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.ob.deg)])
  
  
  for (n in c(1:ncol(DEG.stat))) {
    temp.compair <- colnames(DEG.stat)[n]
    c1=cell.list[[unlist(strsplit(temp.compair,split="_vs_"))[1]]]
    c2=cell.list[[unlist(strsplit(temp.compair,split="_vs_"))[2]]]
    
    data.ob.deg.temp= subset(data.ob.deg,cell=c(c1,c2))
    ident=as.vector(Idents(data.ob.deg.temp))
    names(ident)=names(Idents(data.ob.deg.temp))
    ident[c1]="c1"
    ident[c2]="c2"
    ident=as.factor(ident)
    Idents(data.ob.deg.temp)=ident
    #G1G2.DEG <- suppressMessages(FindMarkers(data.ob.deg.temp,ident.1="c1",ident.2="c2",test.use="MAST",assay="RNA",verbose=F,logfc.threshold=0.1,min.pct = 0.1,pseudocount.use=0.1) %>% tibble::rownames_to_column(var="gene")) %>% tbl_df() %>% left_join(Gene.discrp,by="gene")
    G1G2.DEG <- suppressMessages(FindMarkers(data.ob.deg.temp,ident.1="c1",ident.2="c2",test.use="MAST",assay="RNA",verbose=F,logfc.threshold=0.1,min.pct = 0.1,pseudocount.use=0.1) %>% tibble::rownames_to_column(var="gene")) %>% tbl_df() %>% left_join(Gene.discrp,by="gene")
    
    
    
    # pct.cut off
    pct.cutoff <- 1/10
    G1.sig.up <- subset(G1G2.DEG,G1G2.DEG$avg_log2FC > 0.58 & G1G2.DEG$p_val_adj <0.05 & G1G2.DEG$pct.1 >= pct.cutoff)
    G2.sig.up <- subset(G1G2.DEG,G1G2.DEG$avg_log2FC < -0.58 & G1G2.DEG$p_val_adj <0.05 & G1G2.DEG$pct.2 >= pct.cutoff)
    G1G2.DEG.sig <- rbind(G1.sig.up,G2.sig.up)
    
    temp.up.ID=G1.sig.up$gene
    temp.down.ID=G2.sig.up$gene
    
    DEG.stat["up_regulated",temp.compair]=length(temp.up.ID)
    DEG.stat["down_regulated",temp.compair]=length(temp.down.ID)
    
    print(temp.compair)  
    DEG.results[[temp.compair]] <- list()
    DEG.results[[temp.compair]] [["DEG.all.result"]] <- G1G2.DEG
    DEG.results[[temp.compair]] [["DEG.result"]] <- G1G2.DEG.sig
    DEG.results[[temp.compair]] [["DEG.result.up"]] <-  G1.sig.up 
    DEG.results[[temp.compair]] [["DEG.result.down"]] <-  G2.sig.up 
  }
  save(DEG.results,DEG.stat, file=savefile)
}

#traj.dim.out <- cds.out.DM %>% select(cell:seurat_clusters,Pseudotime,State,cluster_EML_state,traj_Dim1,traj_Dim2) %>% rename(traj_dim1=traj_Dim1,traj_Dim2=traj_Dim2)
#save(DEG.results,traj.dim.out,file="tmp_data/check.traj.DEG.Rdata")





#' check the DEG expression

for(n in colnames(DEG.stat)) {
  print(n)
  temp.compair <- n
  c1=cell.list[[unlist(strsplit(temp.compair,split="_vs_"))[1]]]
  c2=cell.list[[unlist(strsplit(temp.compair,split="_vs_"))[2]]]
  cells.use <- c(c1,c2)
  
  
  up.genes.use <- DEG.results[[temp.compair]]$DEG.result.up$gene
  down.genes.use <- DEG.results[[temp.compair]]$DEG.result.down$gene
  
  colorsR <- cds.out.DM%>% filter(cell %in% c(c1,c2)) %>% select(cell,devTime,cluster_EML,State,cluster_EML_state)  %>%  arrange(State,cluster_EML_state)
  colorsR <- colorsR %>% filter(cell %in% c1) %>% bind_rows(colorsR %>% filter(cell %in% c2)) %>% tibble::column_to_rownames("cell")
  
  if (nrow(DEG.results[[temp.compair]]$DEG.result) != 0) {
    print(paste("Heatmap of gene expression in",temp.compair))
    
    if (length(up.genes.use) >0) {
      if (length(down.genes.use) >0) {
        pheatmap(FunPreheatmapNoLog(lognormExp.mBN,c(up.genes.use,down.genes.use),c(c1,c2)),cluster_rows=F,cluster_cols=F,border_color = "NA",colorRampPalette(c("royalblue3","white","firebrick4"))(50),gaps_col=rep(length(c1),5),gaps_row = rep(length(up.genes.use),5),show_rownames = F,show_colnames = F,main=temp.compair,annotation_col=colorsR )
      }else{
        pheatmap(FunPreheatmapNoLog(lognormExp.mBN,c(up.genes.use,down.genes.use),c(c1,c2)),cluster_rows=F,cluster_cols=F,border_color = "NA",colorRampPalette(c("royalblue3","white","firebrick4"))(50),gaps_col=rep(length(c1),5),show_rownames = F,show_colnames = F,main=temp.compair,annotation_col=colorsR )
      }
    }else{
      pheatmap(FunPreheatmapNoLog(lognormExp.mBN,c(up.genes.use,down.genes.use),c(c1,c2)),cluster_rows=F,cluster_cols=F,border_color = "NA",colorRampPalette(c("royalblue3","white","firebrick4"))(50),gaps_col=rep(length(c1),5),show_rownames = F,show_colnames = F,main=temp.compair,annotation_col=colorsR ) 
    }
  }else{
    print(paste("No DEGs in",temp.compair))
  }
}

#' #### output results
savedir <- paste0("tmp_data/",TD,"/DEG_MAST/")
for (n in c(1:length(colnames(DEG.stat)))) {
  temp.compair <- colnames(DEG.stat)[n]
  temp.out <- DEG.results[[temp.compair]]
  write.table(temp.out$DEG.all.result ,file=paste(savedir,temp.compair,".all.including.nonSig",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  if (nrow(temp.out$DEG.result.up) > 0) {
    write.table(temp.out$DEG.result.up ,file=paste(savedir,temp.compair,".upDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    # write.table(temp.out[["DEG.up.GO.result"]][["out"]],file=paste(savedir,temp.compair,".upDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  }else{
    write.table("None",file=paste(savedir,temp.compair,".upDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    #write.table("None",file=paste(savedir,temp.compair,".upDEG.GO",sep=""),col.names = F,row.names = F,sep="\t",quote=F)
  }
  if (nrow(temp.out$DEG.result.down) > 0) {
    write.table(temp.out$DEG.result.down ,file=paste(savedir,temp.compair,".downDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    #write.table(temp.out[["DEG.down.GO.result"]][["out"]],file=paste(savedir,temp.compair,".downDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  }else{
    write.table("None",file=paste(savedir,temp.compair,".downDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    #write.table("None",file=paste(savedir,temp.compair,".downDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  }
}