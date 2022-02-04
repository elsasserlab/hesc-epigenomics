#' ---
#' title: "Define EZH2I cluster"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---


#' R3.6
# include=FALSE
rm (list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

#' workding directory
DIR="/home/chenzh/My_project/Nerges_EZH2i_ScRNA"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
knitr::opts_chunk$set(echo=FALSE)
setwd(DIR)


#' loading R packages
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(monocle))
suppressMessages(library(ggvenn))
suppressMessages(library(scales))



#' source local function
source("~/PC/SnkM/SgCell.R")
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
#source("src/local.quick.function.R") 

TD="Nov5_2021"


#' loading dataset
data.ob <- readRDS(paste0("tmp_data/",TD,"/multi.IT.mnn.rds"))
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/multi.IT.umap.cord.rds"))%>% mutate(cluster_EML=ifelse(cluster_EML=="ExE_MeLC" & pj=="Nerges_EZH2i_Naive","MeLC",cluster_EML))%>% mutate(rename_EML=ifelse(pj=="Nerges_EZH2i_Naive",cluster_EML,rename_EML)) 
#lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/combineRefH.lognormExp.mBN.rds"))
lognormExp.NP.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.NPonly.rds"))

#cds.out.DM <- readRDS(paste0("tmp_data/",TD,"/naive.trajectory.dm.rds"))
load(paste0("tmp_data/",TD,"/naive.trajectory.Rdata"),verbose = T)
rm(cds,cds.diff_test_res,disp_table)
load(paste0("tmp_data/",TD,"/psd.cutoff.Rdata"),verbose=T)

#' loading traj.DEG results
load(paste0("tmp_data/",TD,"/DEG_branch_MAST.Rdata"),verbose = T)
#' figure setting
#lineage.col.set <-c("Amnion"="#F8766D" ,"AMLC"="#F8766D" ,"EZH2i_Naive_s0"="#EAABA9" ,"EZH2i_Naive_s1"="#F8766D" , "EZH2i_Naive_s2"="#00BADE" , "Epiblast"="#00BADE","Endoderm"="#B385FF","EZH2i_Naive_s3"="#B385FF","TE"= "#64B200", "EZH2i_Naive_s4"= "#64B200","MeLC"= "#EF67EB","Mesoderm"= "#EF67EB", "Prelineage"= "#7B554E", "PriS"="#DB8E00",  "Undef"= "grey50", "ExE_Mes"="darkblue","hES"="#51C0CC","PGC"="#C4423E","hPGCLC"="#C4423E","EZH2i_Primed_D7"="royalblue","EZH2i_Primed_WT"="#00BADE")
lineage.col.set <-c("Amnion"="#E41A1C" ,"AMLC"="#E41A1C" , "ELC"="#4DAF4A" , "Epiblast"="#4DAF4A","Endoderm"="#984EA3","HLC"="#984EA3","TE"= "#377EB8", "TLC"= "#51C0CC","MeLC"= "#FCCDE5","Mesoderm"= "#FCCDE5", "Prelineage"= "#B35806","ICM"="#FFFFB3", "PriS"="#DB8E00",  "Undef"= "grey50", "Unknown"="grey33","ExE_Mes"="#F781BF","hES"="#00FF92","hPGCLC"="#CD8BBC","EarlyBlastocyst"="#FDB863","8C_Morula"="#B35806","Mes"="#FCCDE5","EPI_Amnion"="#FDB863","EZH2i_Primed_D7"="#4DAF4A","EZH2i_Primed_WT"="#4DAF4A","ExE_MeLC"="#F781BF", "gELC"="#D39200","aELC"="#93AA00","TaELC"="#619CFF","MaELC"="#DB72FB")


dt.col.set <- c("EZH2i_Naive_WT"="#7CAE00","EZH2i_Naive_D2"="#F8766D","EZH2i_Naive_D4"="#C77CFF","EZH2i_Naive_D7"="#00BFC4","EZH2i_Primed_WT"="#A6D769","EZH2i_Primed_D7"="#1B984E")

pj.shape.set <- c(CS7=21,D3post=22,SPH2016=24)
pj.shape.set.solid <- c(CS7=16,D3post=15,SPH2016=17)
label.pj <- list()
label.pj$SPH2016="Petropoulos et al., 2016"
label.pj$D3post="Xiang et al., 2020"
label.pj$CS7="Tyser et al., 2021"
label.pj$JPF2019="Zheng et al., 2019"

plot.results <- list() 
plot.results.jpeg <- list()

# text position embryonic dataset
data.EM.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>%group_by(rename_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()  %>% filter(!rename_EML %in% c("3D_ICM","PSA-EPI")) %>% mutate(UMAP_1=ifelse(rename_EML=="PriS",UMAP_1-1,UMAP_1)) %>% mutate(rename_EML=gsub("ExE_Mes","Extraembryonic\nMesoderm",rename_EML))%>% mutate(rename_EML=gsub("PriS","Primitive\nStreak",rename_EML))

# text position of other datasets
data.LC.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("Nerges_EZH2i_Naive","Nerges_EZH2i_Primed","JPF2019"))  %>%group_by(rename_EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup() 




plot.results$em <- ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("D3post","CS7","SPH2016")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>% mutate(od=factor(rename_EML,c("Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","PGC","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=rename_EML,shape=pj),size=2,color="grey")+ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2021"))+NoAxes()+NoLegend()

umap.xlim=ggplot_build(plot.results$em)$layout$panel_scales_x[[1]]$range$range
umap.ylim=ggplot_build(plot.results$em)$layout$panel_scales_y[[1]]$range$range

plot.results.jpeg$umap.em  <-ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("D3post","CS7","SPH2016")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>% mutate(od=factor(rename_EML,c("Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","PGC","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=rename_EML,shape=pj),size=2,color="grey")+ggtitle("Embryonic dataset")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2021"))+NoAxes()+NoLegend()

plot.results$em_legend <- ggplot()+ geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016")) %>% filter(!rename_EML %in% c("PSA-EPI","3D_ICM")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML,shape=pj),size=2)+theme_classic()+scale_color_manual("",values=lineage.col.set)+scale_shape_manual("",values = pj.shape.set.solid,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2020",'CS7'="Tyser et al., 2020")) + theme(legend.text = element_text( size = 6),legend.title = element_text( size = 1))
plot.results$em_legend <- plot.results$em_legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()


plot.results$JPF <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.25)+ggtitle("Zheng et al., 2019")+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE","PGC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="JPF2019"),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$JPF <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.25)+ggtitle("Zheng et al., 2019")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)

plot.results$em_JPF_legend <- ggplot() +geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=2)+ggtitle("Zheng et al., 2019") +theme_classic()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual("",values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
plot.results$em_JPF_legend <- plot.results$em_JPF_legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()


#plot.results$Nerges_EZH2i_Naive <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("Nerges_EZH2i_Naive")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("Nerges_EZH2i_Naive")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.75)+ggtitle("Nerges_EZH2i_Naive")+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="Nerges_EZH2i_Naive"),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
# 
#plot.results$Nerges_EZH2i_Primed <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("Nerges_EZH2i_Primed")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("Nerges_EZH2i_Primed")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.75)+ggtitle("Nerges_EZH2i_Primed")+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE","PGC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="Nerges_EZH2i_Primed"),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)


plot.results$Nerges_EZH2i_Naive_devTime <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("Nerges_EZH2i_Naive")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("Nerges_EZH2i_Naive")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=devTime),size=0.75)+ggtitle("Nerges_EZH2i_Naive")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=dt.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
plot.results$Nerges_EZH2i_Primed_devTime <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("Nerges_EZH2i_Primed")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("Nerges_EZH2i_Primed")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=devTime),size=0.75)+ggtitle("Nerges_EZH2i_Primed")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=dt.col.set)+xlim(umap.xlim)+ylim(umap.ylim)

plot.results$devTime_legend <-ggplot()+geom_point( data.ob.umap %>% filter(pj %in% c("Nerges_EZH2i_Primed","Nerges_EZH2i_Naive")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=devTime),size=2)+ggtitle("")+ theme_classic()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual("",values=dt.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
plot.results$devTime_legend <- plot.results$devTime_legend  %>% ggpubr::get_legend() %>% ggpubr::as_ggplot()



plot.results$Nerges_EZH2i_Naive <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("Nerges_EZH2i_Naive")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("Nerges_EZH2i_Naive")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.75)+ggtitle("Nerges_EZH2i_Naive")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
# 
plot.results$Nerges_EZH2i_Primed <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("Nerges_EZH2i_Primed")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("Nerges_EZH2i_Primed")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.75)+ggtitle("Nerges_EZH2i_Primed")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)




for (n in c("EZH2i_Naive_WT","EZH2i_Naive_D2","EZH2i_Naive_D4","EZH2i_Naive_D7","EZH2i_Primed_WT","EZH2i_Primed_D7")) {
  plot.results[[n]] <-  ggplot()+ geom_point(data.ob.umap %>% filter(devTime != n) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(devTime==n) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.75)+ggtitle(gsub("EZH2i_","",n))+ theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)#+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")
  plot.results[[paste(n,"_high",sep="")]] <-  ggplot()+ geom_point(data.ob.umap %>% filter(devTime != n) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(devTime==n) ,mapping=aes(x=UMAP_1,y=UMAP_2),size=0.25,col="red")+ggtitle(gsub("EZH2i_","",n))+ theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)#+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")
}

plot.results$Naive_legend <-ggplot()+geom_point( data.ob.umap %>% filter(devTime=="EZH2i_Naive_D7") ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=2)+ggtitle(gsub("EZH2i_","",n))+ theme_classic()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual("",values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
plot.results$Naive_legend <- plot.results$Naive_legend %>% ggpubr::get_legend() %>% ggpubr::as_ggplot()


data.cluster.text.pos <- data.ob.umap %>%group_by(seurat_clusters) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()
plot.results$umap.cluster <- ggplot()+ geom_point(data.ob.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters),size=0.5)+geom_text(data.cluster.text.pos  ,mapping=aes(x=UMAP_1,y=UMAP_2,label=seurat_clusters),size=4,color="Black")+theme_classic()+NoAxes()+NoLegend()+ggtitle("Seurat clusters")+theme(plot.title = element_text(hjust=0.5))#+scale_color_manual(values=cluster.col.set)
plot.results.jpeg$umap.cluster <- ggplot()+ geom_point(data.ob.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters),size=0.5)+theme_classic()+NoAxes()+NoLegend()+ggtitle("Seurat clusters")+theme(plot.title = element_text(hjust=0.5))#



#' output the number of different cells
data.ob.umap  %>% filter(pj=="Nerges_EZH2i_Naive")%>% group_by(pj,devTime,cluster_EML) %>% summarise(nCell=n()) %>% ungroup() %>% arrange(devTime) %>% select(-pj) %>% write.table( paste0("tmp_data/",TD,"/cell.number.tsv"),sep="\t",quote=F,col.names = T,row.names = F)

temp <- data.ob.umap  %>% filter(pj=="Nerges_EZH2i_Naive") %>% mutate(cla=cluster_EML) %>% group_by(pj,devTime,cla) %>% summarise(nCell=n()) %>% group_by(pj,devTime) %>% mutate(prop=nCell/sum(nCell)) %>% mutate(devTime=gsub("EZH2i_Naive_","",devTime)) %>% mutate(devTime=ifelse(devTime=="WT","D0",devTime)) %>% mutate(devTime=factor(devTime,c("D0","D2","D4","D7"),ordered = T)) %>% ungroup() %>% arrange(devTime)

plot.results$Naive_nCell <- temp %>% ggplot+geom_bar(mapping=aes(x=devTime,fill=cla,y=prop*100),stat="identity",position="dodge")+ theme_classic() +xlab("")+ylab("% cells")+ggtitle("Nerges_EZH2i_Naive")+FunTitle()+scale_fill_manual("",values=lineage.col.set[c("ELC","HLC","TLC","AMLC","MeLC","Undef")])

#‘ check some typical marker gene expression
data.temp <- data.ob %>% subset(cell=(data.ob.umap %>% filter(pj=="Nerges_EZH2i_Naive") %>% filter(!(rename_EML=="Undef") & (pj=="Nerges_EZH2i_Naive")) %>% pull(cell)))
Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"rename_EML"],c("ELC","HLC","TLC", "MeLC", "AMLC"),ordered=T)

#sel.features <- c("DPPA5","POU5F1","NANOG","GATA4","GATA6","PDGFRA","GATA2","GATA3","DLX3","LIX1","TMEM88","PMP22","LUM","POSTN","COL3A1","KRT19","RBPMS","GABRP")
sel.features <- c("DPPA5","POU5F1","NANOG","GATA4","GATA6","PDGFRA","GATA2","GATA3","DLX3","LIX1","TMEM88","PMP22","KRT19","RBPMS","GABRP")
plot.results$Naive_mk <- DotPlot(data.temp,features = sel.features)+ theme(axis.text.x=element_text(angle = 90))+xlab("Genes")+ylab("")

for ( g in c(sel.features)) {
  plot.results[[paste("mk",g,sep="_")]] <- VlnPlot(data.temp,g)+scale_fill_manual("",values=lineage.col.set[c("ELC","HLC","TLC","AMLC","MeLC","Undef")])+NoLegend()+xlab("")
}

#' plot trajectory
plot.results$traj.cluster_EML <- plot_cell_trajectory(cds.out, color_by = "cluster_EML",show_branch_points=F,cell_size=0.75)+NoAxes()+NoLegend()+scale_color_manual(values=lineage.col.set[c("Undef","ELC","HLC","TLC","MeLC","AMLC","Undef")])
plot.results$traj.cluster_EML.ld <- plot_cell_trajectory(cds.out, color_by = "cluster_EML",show_branch_points=F)+NoAxes()scale_color_manual(values=lineage.col.set[c("Undef","ELC","HLC","TLC","MeLC","AMLC","Undef")])

plot.results$traj.state <- plot_cell_trajectory(cds.out, show_branch_points=F,cell_size=0.75)+NoAxes()+NoLegend()
plot.results$traj.state.ld <- plot_cell_trajectory(cds.out, show_branch_points=F)+NoAxes()

plot.results$traj.cluster_EML_state <-plot_cell_trajectory(cds.out, color_by = "cluster_EML_state",show_branch_points=F,cell_size=0.75)+NoAxes()+NoLegend()+scale_color_manual(values=lineage.col.set[c("Undef","aELC","gELC","MaELC","TaELC","HLC","TLC","MeLC","AMLC","Undef")])
plot.results$traj.cluster_EML_state.ld <-plot_cell_trajectory(cds.out, color_by = "cluster_EML_state",show_branch_points=F)+NoAxes()+scale_color_manual(values=lineage.col.set[c("Undef","aELC","gELC","MaELC","TaELC","HLC","TLC","MeLC","AMLC","Undef")])

plot.results$traj.psd <-plot_cell_trajectory(cds.out, color_by = "Pseudotime",show_branch_points=F,cell_size=0.75)+NoAxes()+NoLegend()+ scale_colour_viridis_c()
plot.results$traj.psd.ld <-plot_cell_trajectory(cds.out, color_by = "Pseudotime",show_branch_points=F)+NoAxes() + scale_colour_viridis_c()

plot.results$psd.dis <- cds.out.DM %>% mutate(od=factor(devTime,c("EZH2i_Naive_WT","EZH2i_Naive_D2","EZH2i_Naive_D4","EZH2i_Naive_D7"),ordered = T)) %>% arrange(od) %>% mutate()%>% ggplot(mapping=aes(x=Pseudotime,y = od,fill=devTime))+ggridges::geom_density_ridges()+theme_classic()+ggtitle("EZH2i_Naive cells")+FunTitle()+ylab("")+NoLegend()+scale_fill_manual(values=dt.col.set)+
geom_vline(xintercept = psd.cutoff,linetype="dashed")+geom_vline(xintercept = psd.branch,linetype="dashed")
branch.layer <- plot.results$traj.state$layer[[1]]

# highlight cells on trajectory
temp.plot <- list()
for (n in unique(cds.out.DM$devTime)) {
  temp.input <- cds.out.DM %>% mutate(sel=ifelse(devTime==n , "sel","no")) %>% arrange(sel)
  temp.plot[[paste("traj.hl",n,sep=".")]] <-  ggplot()+branch.layer+geom_point(data=(temp.input %>% filter(sel=="no")),mapping=aes(x=traj_Dim1,y=traj_Dim2,col=selected),size=0.5,col="grey")+geom_point(data=(temp.input %>% filter(sel=="sel")),mapping=aes(x=traj_Dim1,y=traj_Dim2,col=selected),size=0.5,col="red")+theme_classic()+NoAxes()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5))
}
cowplot::plot_grid(plotlist=temp.plot[c(3,1,2,4)])
plot.results$traj.hl.tp <- temp.plot[c(3,1,2,4)]

temp.plot <- list()
for (n in unique(cds.out.DM$cluster_EML)) {
  temp.input <- cds.out.DM %>% mutate(sel=ifelse(cluster_EML==n , "sel","no")) %>% arrange(sel)
  temp.plot[[paste("traj.hl",n,sep=".")]] <-  ggplot()+branch.layer+geom_point(data=(temp.input %>% filter(sel=="no")),mapping=aes(x=traj_Dim1,y=traj_Dim2,col=selected),size=0.5,col="grey")+geom_point(data=(temp.input %>% filter(sel=="sel")),mapping=aes(x=traj_Dim1,y=traj_Dim2,col=selected),size=0.5,col="red")+theme_classic()+NoAxes()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5))
}
cowplot::plot_grid(plotlist=temp.plot[c(1,4,2,3,6,5)])
plot.results$traj.hl.cluster <- temp.plot[c(1,4,2,3,6,5)]

#' check the DEG overlap
temp.list <- list(
  D7_TaELC_vs_D7_MaELC = DEG.results$D7_TaELC_vs_D7_MaELC$DEG.result %>% mutate(UD=ifelse(avg_log2FC> 0,"Up","Down")) %>% mutate(ID=paste(gene,UD)) %>% pull(ID) ,
  D7_TLC_vs_D7_MeLC = DEG.results$D7_Naive_TLC_vs_D7_Naive_MeLC$DEG.result %>% mutate(UD=ifelse(avg_log2FC> 0,"Up","Down")) %>% mutate(ID=paste(gene,UD)) %>% pull(ID)
)
plot.results$TvM <- ggvenn(temp.list,set_name_size=4)


#' check the Vlnplot of Traj DEG
temp.sel.gene <- c("HAND1","UTF1","LEFTY2","TFAP2C","NODAL")
temp.M <- cds.out.DM %>% filter(cluster_EML %in% c("ELC","MeLC","TLC")) %>% filter(devTime=="EZH2i_Naive_D7")
temp.sel.exp <- lognormExp.NP.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(cluster_EML_state,c("gELC","aELC","MaELC","TaELC","MeLC","TLC"),ordered = T)) %>% arrange(od)  

temp.sel.gene <- c("HAND1","UTF1","LEFTY2","TFAP2C","NODAL")
for (g in temp.sel.gene ) {
  plot.results[[paste0("traj.DEG.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=cluster_EML_state),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")+scale_fill_manual(values=lineage.col.set[c("gELC","aELC","MaELC","TaELC","TLC","MeLC")])
  
  plot.results.jpeg[[paste0("traj.DEG.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=cluster_EML_state),color="black",scale = "width")+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")+scale_fill_manual(values=lineage.col.set[c("gELC","aELC","MaELC","TaELC","TLC","MeLC")])
  
}

#' check the gene expression along the traj
temp.gene <- c("NANOG","POU5F1","SOX2","GATA3","TBXT","HAND1","TFAP2C","LEFTY2","NODAL")
temp.input <- lognormExp.NP.mBN[temp.gene,] %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% gather(cell,logExp,-gene) %>% inner_join(cds.out.DM,by="cell") 

temp.plot <- list()
for (g in temp.gene) {
  temp.plot[[g]] <- ggplot()+geom_smooth(data=(temp.input %>% filter(gene==g) %>% filter(State %in% c(1,2))),mapping=aes(x=Pseudotime,y=logExp),col=lineage.col.set["TLC"],method="loess")+geom_smooth(data=(temp.input %>% filter(gene==g) %>% filter(State %in% c(1,3))),mapping=aes(x=Pseudotime,y=logExp),col=lineage.col.set["MeLC"],method="loess")+theme_classic()+ggtitle(g)+FunTitle()+NoLegend()+geom_vline(xintercept = psd.cutoff ,linetype="dashed")+geom_vline(xintercept = psd.branch ,linetype="dashed")+xlim(0,18)+scale_y_continuous(limit=c(0,NA),oob=squish)
}
cowplot::plot_grid(plotlist=temp.plot)
pdf("tmp_data/figure.traj.exp.pdf",9,9)
cowplot::plot_grid(plotlist=temp.plot)
dev.off()

pdf("tmp_data/figure.UMAP.related.pdf")
plot.results$em
plot.results.jpeg$umap.em
plot.results$em_legend
plot.results$JPF
plot.results.jpeg$JPF
plot.results$em_JPF_legend 

plot.results$umap.cluster
plot.results.jpeg$umap.cluster

plot.results$Nerges_EZH2i_Naive_devTime
plot.results$Nerges_EZH2i_Primed_devTime
plot.results$devTime_legend
plot.results$Nerges_EZH2i_Naive
cowplot::plot_grid(
  plot.results$EZH2i_Naive_WT,
  plot.results$EZH2i_Naive_D2,
  plot.results$EZH2i_Naive_D4,
  plot.results$EZH2i_Naive_D7
)
plot.results$Naive_legend 
cowplot::plot_grid(
  plot.results$EZH2i_Naive_WT_high,
  plot.results$EZH2i_Naive_D2_high,
  plot.results$EZH2i_Naive_D4_high,
  plot.results$EZH2i_Naive_D7_high
)
plot.results$Nerges_EZH2i_Primed
cowplot::plot_grid(
  plot.results$EZH2i_Primed_WT,
  plot.results$EZH2i_Primed_D7,
  nrow=2,
  ncol=2
)
cowplot::plot_grid(
  plot.results$EZH2i_Primed_WT_high,
  plot.results$EZH2i_Primed_D7_high,
  nrow=2,
  ncol=2
)
dev.off()

pdf("tmp_data/figure.nCell.bar.pdf",5.5,4)
plot.results$Naive_nCell
dev.off()
pdf("tmp_data/DotPlot.mk.EZH2i_Naive.bar.pdf",7.5,4)
plot.results$Naive_mk 
dev.off()

pdf("tmp_data/VlnPlot.mk.EZH2i_Naive.bar.pdf",12,12)
cowplot::plot_grid(
  plotlist=plot.results[paste("mk",sel.features,sep="_")]
)
dev.off()

pdf("tmp_data/figure.traj.pdf",12,4)
cowplot::plot_grid(
  plot.results$traj.state,
  plot.results$traj.cluster_EML_state,
  plot.results$traj.psd,
  plot.results$traj.cluster_EML,nrow=1
)
cowplot::plot_grid(
  plot.results$traj.state.ld,
  plot.results$traj.cluster_EML_state.ld,
  plot.results$traj.psd.ld,
  plot.results$traj.cluster_EML.ld
)
dev.off()

pdf("tmp_data/figure.traj.hl1.pdf",12,4)
print(
  cowplot::plot_grid(plotlist=plot.results$traj.hl.tp,nrow=1)
)
dev.off()
pdf("tmp_data/figure.traj.hl2.pdf",9,8)
print(
  cowplot::plot_grid(plotlist=plot.results$traj.hl.cluster,nrow=2)
)
dev.off()

pdf("tmp_data/DEG.TvM.ov.pdf",4,4)
plot.results$TvM
dev.off()

pdf("tmp_data/Psd.hist.plot.pdf",6,4)
plot.results$psd.dis
dev.off()

pdf("tmp_data/figure.traj.DEG.pdf",9,8)
print(
  cowplot::plot_grid(plotlist=plot.results[paste0("traj.DEG.Vln.",temp.sel.gene)])
)
print(
  cowplot::plot_grid(plotlist=plot.results.jpeg[paste0("traj.DEG.Vln.",temp.sel.gene)])
)
dev.off()







