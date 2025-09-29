setwd("/work/wjx/scRNA-seq/AIP/BCR/filter_contig_annotations/")
library(scRepertoire)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(patchwork)
library(ggsci)
library(ggalluvial)


PatientList <- grep(pattern = "filtered_contig_annotations",list.files("/work/wjx/scRNA-seq/AIP/BCR/filter_contig_annotations"),value=T)
WholeTissue_bcr_list <- lapply(PatientList,function(x){
  cellclone <- read.csv(paste("/work/wjx/scRNA-seq/AIP/BCR/filter_contig_annotations/",x,sep = ""))
  cellclone <-cellclone[,c("barcode","contig_id","length","chain","v_gene", "d_gene","j_gene","c_gene","cdr3","cdr3_nt","reads","umis","raw_clonotype_id")]
})
names(WholeTissue_bcr_list)=gsub("_filter.*","",PatientList)
WholeTissue_bcr=do.call(rbind, WholeTissue_bcr_list)


WholeTissue_bcr$orig.ident=gsub("[.].*","",rownames(WholeTissue_bcr))
WholeTissue_bcr$Patient=gsub("P1.*","P1",WholeTissue_bcr$orig.ident) 
WholeTissue_bcr$Status=gsub("P.*","AIP",WholeTissue_bcr$orig.ident)
WholeTissue_bcr$Status=gsub("N.*","Normal",WholeTissue_bcr$Status)


ClusterFreq <- WholeTissue_bcr %>% dplyr::filter(chain %in% c("IGH")) %>% dplyr::select("c_gene","Status") %>% table %>%
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterFreq <- dplyr::filter(ClusterFreq,grepl("IGH",seurat_clusters)) 
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters <- factor(ClusterPer$seurat_clusters,levels = c("IGHM","IGHD","IGHE","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4"))
pdf("SF3a.pdf",width = 10,height = 5)
ggplot(data = na.omit(ClusterPer), aes(x = stim, fill=seurat_clusters,y=Per)) + 
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=c("#6495ED","#0000CD","#9400D3","#90EE90","#32CD32","#FFD700","#FFA500","#FF69B4","#DC143C"))+
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 17),axis.text.x=element_text(color = "black",size = 17),
        axis.title.x = element_text(size = 20),axis.title.y  = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=13))
dev.off()

ClusterFreq <- WholeTissue_bcr %>% dplyr::filter(chain %in% c("IGH")) %>% dplyr::select("c_gene","Patient") %>% table %>%
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterFreq <- dplyr::filter(ClusterFreq,grepl("IGH",seurat_clusters)) 
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters <- factor(ClusterPer$seurat_clusters,levels = c("IGHM","IGHD","IGHE","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4"))
ClusterPer$stim <- factor(ClusterPer$stim,levels = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","N1","N2","N3","N4","N5","N6","N7","N8"))
pdf("SF3b.pdf",width = 10,height = 10)
ggplot(data = na.omit(ClusterPer), aes(x = stim, fill=seurat_clusters,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=c("#6495ED","#0000CD","#9400D3","#90EE90","#32CD32","#FFD700","#FFA500","#FF69B4","#DC143C"))+
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 17),axis.text.x=element_text(color = "black",size = 17),
        axis.title.x = element_text(size = 20),axis.title.y  = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=13))
dev.off()





DefineB <- readr::read_rds("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")

PatientList <- list.files("/work/wjx/scRNA-seq/AIP/BCR/filter_contig_annotations")

Mtcr_list <- lapply(PatientList,function(x){
  Xdata <- read.csv(file = paste("/work/wjx/scRNA-seq/AIP/BCR/filter_contig_annotations/",x,sep = ""), stringsAsFactors = F)
  return(Xdata)
})
names(Mtcr_list)=gsub("_filtered_contig_annotations.csv","",PatientList )
head(Mtcr_list[[1]])

barcode <- rownames(DefineB@meta.data)


combined_Mtcr <- combineBCR(Mtcr_list,samples = names(Mtcr_list),ID= names(Mtcr_list),removeNA=T,call.related.clones=T,threshold = 0.85)
names(combined_Mtcr) <- names(Mtcr_list)
head(combined_Mtcr[[1]])
combined_Mtcr <- lapply(names(combined_Mtcr),function(Name){
  x <- combined_Mtcr[[Name]]
  rx <- Mtcr_list[[Name]]
  x$barcode <- gsub("N[1-8]_N[1-8]_|P[1-9]_P[1-9]_|P1_1_P1_1_|P1_2_P1_2_","",x$barcode)
  x$raw_clonotype_id  <- rx$raw_clonotype_id[match(x$barcode,rx$barcode)]
  x$barcode <- paste(Name,"AIP_Inflam",x$barcode,sep="_")
  if(Name %in% c("N1","N2","N3","N4","N5","N6","N7","N8")){
    x$barcode <- gsub("Inflam","Normal",x$barcode)
  }
  rownames(x) <- x$barcode
  return(x)
})
names(combined_Mtcr) <- names(Mtcr_list)
head(combined_Mtcr[[1]])
combined_Mtcr <- addVariable(combined_Mtcr, variable.name = "patient", variables = c("N1","N2","N3","N4","N5","N6","N7","N8","P1","P1","P2","P3","P4","P5","P6","P7","P8","P9"))
combined_Mtcr <- addVariable(combined_Mtcr, variable.name = "group", variables = c(rep("Normal",8),rep("AIP",10)))
readr::write_rds(combined_Mtcr,path="combined_Mtcr_20240715_callrelated.rds.gz",compress="gz")
combined_Mtcr <- readRDS("/work/wjx/scRNA-seq/AIP/BCR/combined_Mtcr_20240715_callrelated.rds.gz")


data <- clonalQuant(combined_Mtcr, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE,
            group.by = "patient",
            exportTable = T,
            palette = "viridis")
data$group <- gsub("N.*","Normal",gsub("P.*","AIP",data$patient))
pdf('SF3d.pdf',width = 4,height = 3.5)
ggplot(data = data,aes(x=group,color=group,y=scaled))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_jitter(width=0.2,size=2.5)+
  ggtitle("Unique Clones")+
  scale_x_discrete(limits = c('Normal',"AIP"))+
  scale_y_continuous(limits = c(NA,max(data$scaled)*1.05))+
  scale_color_manual(limits=c("AIP","Normal"),values=c("#CB181D","#1f78b4"))+
  ylab('Percentage (%)')+
  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(color = 'black'),panel.background = element_blank(),panel.grid = element_blank(),
     legend.title = element_blank(),axis.text.y = element_text(color = 'black'),strip.background = element_blank(),
     strip.text = element_text(size = 12),axis.title.x = element_blank(),axis.line = element_line(colour = 'black'),legend.position = "none")+
  ggpubr::stat_compare_means(comparisons = list(c('Normal','AIP')),paired = F,method = 'wilcox.test')
dev.off()



library(ggh4x)
data <- clonalDiversity(combined_Mtcr, 
                cloneCall = "strict",
                group.by="patient",
                metrics = c("shannon", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE"),
                exportTable = T)
data$group <- gsub("N.*","Normal",gsub("P.*","AIP",data$patient))
data <- data[,-8]
data_new <- gather(data,key="diversity",value="y",-"patient",-"group")

pdf('SF3c.pdf',width = 10,height = 7)
ggplot(data = data_new,aes(x=group,color=group,y=y))+
  geom_boxplot(outlier.shape = NA,size=1)+
  geom_jitter(width=0.2,size=2.5)+
  facet_wrap(~diversity,nrow =2,scales = 'free')+
  scale_x_discrete(limits = c('Normal',"AIP"))+
  scale_color_manual(limits=c("AIP","Normal"),values=c("#CB181D","#1f78b4"))+
  ylab('Percentage (%)')+
  theme(plot.title = element_text(hjust = 0.5),axis.text = element_text(color = 'black'),panel.background = element_blank(),panel.grid = element_blank(),
        legend.title = element_blank(),axis.text.y = element_text(color = 'black'),strip.background = element_blank(),
        strip.text = element_text(size = 12),axis.title.x = element_blank(),axis.line = element_line(colour = 'black'),legend.position = "none")+
  ggpubr::stat_compare_means(comparisons = list(c('Normal','AIP')),paired = F,method = 'wilcox.test')
  
dev.off()



all_clone <- combined_Mtcr%>%bind_rows()
table(DefineB$celltype1)


DefineT <- combineExpression(combined_Mtcr, 
                             DefineB, 
                             filterNA = T,
                                   cloneCall="strict", 
                             group.by = "sample", 
                                   proportion = F, 
                                   cloneSize=c(Small=1, Medium=5, Large=20))

tiff("SF3e.tiff",width = 1000,height = 700,res=150)
Seurat::DimPlot(DefineT, group.by = "cloneSize") +
  scale_color_manual(values=c("#DC050C","#FF7F00","#7BAFDE"))
dev.off()


library(ggraph)
.cluster_cols <- c("#B2DF8A","#B17BA6","#FF7F00", "#DC050C","#E7298A","#1965B0","#55A1B1")
pdf("Figure2h.pdf",width = 8,height = 6)
clonalNetwork(DefineT, 
              reduction = "umap", 
              group.by = "celltype1",
              filter.clones = NULL,
              filter.identity = c("Plasma G"),
              filter.graph = T,
              cloneCall = "strict")+
  scale_color_manual(values=.cluster_cols)
dev.off()


pdf("SF3f.pdf",width = 8,height = 6)
clonalOccupy(DefineT, 
               x.axis = "celltype1", 
             proportion = T, 
             label = F)+
  scale_fill_manual(values=c("#DC050C","#FF7F00","#7BAFDE"))
dev.off()


.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
CTstrict_data <- spread(data.frame(table(DefineT$celltype1,DefineT$CTstrict)),Var1,Freq)
overlapclones <- subset(CTstrict_data,CTstrict_data$`ABC(IgD-)`!=0&CTstrict_data$`Plasma G`!=0)
pdf("Figure2i.pdf",width = 8,height = 6)
clonalCompare(DefineT, 
              cloneCall="strict", 
              chain = "both",
              samples = c("ABC(IgD-)", "Plasma G"), 
              clones = overlapclones$Var2,
              top.clones = 5,
              highlight.clones = NULL,
              relabel.clones = F,
              group.by = "celltype1",
              graph = "alluvial",
              exportTable = F,
              palette = "RdYlGn")+
  scale_fill_manual(values = .cluster_cols)
dev.off()



pdf("SF3g.pdf",width = 8,height = 6)
clonalOverlap(DefineT, 
              cloneCall = "aa", 
              method = "cosine",
              group.by="celltype1",
              palette="Viridis")+
  scale_fill_gradientn(colors =  colorRampPalette(c(RColorBrewer::brewer.pal(8,'Reds')),space='rgb')(10),na.value = "white")
dev.off()

