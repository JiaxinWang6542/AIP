setwd("/work/wjx/scRNA-seq/AIP/figure/figure2")
library(Seurat)
library(clustree)
library(ggplot2)
library(magrittr)
library(harmony)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ComplexHeatmap)
library(ggsci)
library(circlize)
library(RColorBrewer)
options(stringsAsFactors = F)
set.seed(42)

WholeTissueData_ident=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure1/WholeTissueData_ident.rds.gz")
table(WholeTissueData_ident$celltype)
WholeTissueData_B_raw <- subset(WholeTissueData_ident,celltype %in% c("B","Plasma"))
table(WholeTissueData_B_raw@meta.data$celltype)


WholeTissueData_B_raw  <- NormalizeData(WholeTissueData_B_raw) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_B_raw)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_B_raw <- ScaleData(WholeTissueData_B_raw,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissueData_B_raw,ndims = 50)
WholeTissueData_B_raw <- FindNeighbors(WholeTissueData_B_raw,dims = 1:20)
WholeTissueData_B_raw  <- RunUMAP(object = WholeTissueData_B_raw, dims = 1:20) 
WholeTissueData_B_raw  <- FindClusters(WholeTissueData_B_raw, resolution = c(seq(0.1,1.8,0.1)), algorithm = 1)

clustree(WholeTissueData_B_raw@meta.data,prefix = "RNA_snn_res.")

Idents(WholeTissueData_B_raw) <- 'RNA_snn_res.0.7'
wholeTissueData_B_raw_0.7_markers <- FindAllMarkers(WholeTissueData_B_raw,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#remove doublet
WholeTissueData_B <- subset(WholeTissueData_B_raw,RNA_snn_res.0.7!=6 & RNA_snn_res.0.7!=10& RNA_snn_res.0.7!=11 & RNA_snn_res.0.7!=14 & RNA_snn_res.0.7!=16 & RNA_snn_res.0.7!=17& RNA_snn_res.0.7!=15)

WholeTissueData_B  <- NormalizeData(WholeTissueData_B) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_B)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_B <- ScaleData(WholeTissueData_B,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_B,ndims = 50)

WholeTissueData_B <- FindNeighbors(WholeTissueData_B,dims = 1:17)
WholeTissueData_B  <- RunUMAP(object = WholeTissueData_B, dims = 1:17) 
WholeTissueData_B  <- FindClusters(WholeTissueData_B, resolution = c(seq(0.1,1.8,0.1)), algorithm = 1)

clustree(WholeTissueData_B@meta.data,prefix = "RNA_snn_res.")

table(WholeTissueData_B$RNA_snn_res.0.6,WholeTissueData_B$orig.ident)

levels(Idents(WholeTissueData_B))
Idents(WholeTissueData_B) <- 'RNA_snn_res.0.8'
WholeTissueData_B_markers_0.8 <- FindAllMarkers(WholeTissueData_B, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


WholeTissueData_B_ident <- RenameIdents(WholeTissueData_B, "0" = "Naive B","1" = "Naive B","7" = "Naive B","2" = "Memory B","4" = "Memory B","8" = "Memory B","9"="ABC(IgD+)","10"="ABC(IgD-)","13"="GC B","3"="Plasma G","5"="Plasma G","6"="Plasma G","12"="Plasma G","11"="Plasma A")



WholeTissueData_B_ident$celltype1=WholeTissueData_B_ident@active.ident
table(WholeTissueData_B_ident$celltype1,WholeTissueData_B_ident$orig.ident)
table(WholeTissueData_B_ident$Status)

.cluster_cols <- c("#B2DF8A","#1965B0","#FB8072", "#DC050C","#FF7F00","#882E72","#B17BA6")
  

readr::write_rds(WholeTissueData_B_ident,path="./result/WholeTissueData_B_ident.rds.gz",compress="gz")

B_markers <- FindAllMarkers(WholeTissueData_B_ident,only.pos = TRUE, min.pct = 0.1,logfc.threshold = 0.25)



WholeTissueData_B_ident$Status <- factor(WholeTissueData_B_ident$Status,levels = c("Normal","AIP"))


tiff('Figure2A.tiff',width = 1400,height = 700,res=150)
DimPlot(WholeTissueData_B_ident,reduction = 'umap', label = F,repel = T,cols = .cluster_cols,group.by = "celltype1",split.by = "Status",pt.size = 0.5)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

WholeTissueData_B_ident$orig.ident <- factor(WholeTissueData_B_ident$orig.ident,levels = c("N1","N2","N3","N4","N5","N6","N7","N8","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))

tiff('SF2A.tiff',width = 1120,height = 910,res = 150)
DimPlot(WholeTissueData_B_ident,reduction = 'umap', label = F,repel = T, group.by = "orig.ident")+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


tiff("SF2B.tiff",width = 1500,height = 1500,res =200)
FeaturePlot(WholeTissueData_B_ident,c("MS4A1","MZB1","TCL1A","AIM2","ITGAX","LMO2","IGHG1","IGHG4","IGHA2"),cols = rev(colorRampPalette(c(RColorBrewer::brewer.pal(10,'Spectral')[2:10]),space='rgb')(20)),ncol = 3)&
  theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line =element_blank())
dev.off()

Idents(WholeTissueData_B_ident) <- "celltype1"





PlotGenes <- c("CD200","COL19A1","TSPAN13","BTLA","CALHM6","TCL1A","FCER2","IL4R","IGHM","BACH2",
               "LINC01781","MARCKS","AIM2","CD24","CD70","LTB","GPR183","COCH","TNFRSF13B","CD27",
               "ITGAX","TBX21","LILRB2","SLC25A29","TENT5A","MPP6","RIN3","LILRA4","PPP1R14A","IGHD",
               "ITGAX","TBX21","AICDA","DHRS9","GPR34","SDR16C5","CCR1","ZBED2","PTMS","TSPO",
               "RGS13","BCL6","SERPINA9","AC023590.1","MEF2B","ELL3","MYO1E","LCK","CD83","LMO2",
               "MZB1","SDC1","IGHG4","SPAG4","CST3","DPEP1","IGKV1-12","IGKV1-16","IGHV5-51","IGKV1-27",
               "MZB1","SDC1","IGHA2","LINC02362","ZBP1","ABCB9","JCHAIN","PRDM1","GPRC5D","LINC02384")

Sub_Data <- WholeTissueData_B_ident

Sub_TF_Matrix <- Sub_Data@assays$RNA@data[PlotGenes,] %>% as.matrix()

Sub_TF_Matrix <- Sub_TF_Matrix[apply(Sub_TF_Matrix,1,sum) > 0,]

Sub_ClusterID <- lapply(split(Sub_Data@meta.data,list(Sub_Data$celltype1)), function(x)rownames(x))

Sub_TF_ClusterMean <- t(apply(Sub_TF_Matrix,1,function(x){
  lapply(Sub_ClusterID,function(ID){
    mean(x[ID])
  }) %>% unlist
}))

Sub_TF_ClusterMean <- Sub_TF_ClusterMean[apply(Sub_TF_ClusterMean,1,sum)>0,]
Sub_TF_ClusterMean[,1:ncol(Sub_TF_ClusterMean)] <- t(apply(Sub_TF_ClusterMean,1,scale))

Sub_barcolor <- .cluster_cols
names(Sub_barcolor) <- unique(Sub_Data$celltype1)
Sub_SampleBar <-  unique(Sub_Data$celltype1)
names(Sub_SampleBar) <-  unique(Sub_Data$celltype1)


library(ComplexHeatmap)
Sub_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Sub_SampleBar,
                                                col = list(bar = Sub_barcolor),show_legend = F,show_annotation_name = F)

pdf("Figure2B.pdf",width = 8,height = 10)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),
                        show_row_names = F,show_column_names = T  ,column_names_rot = 45,top_annotation = Sub_ha_bar)
dev.off()




EachSampleClusterDisB <- lapply(unique(WholeTissueData_B_ident@meta.data$orig.ident),function(i){
  seur <- subset(WholeTissueData_B_ident,cells=rownames(WholeTissueData_B_ident@meta.data[WholeTissueData_B_ident@meta.data$orig.ident %in% i,]))
  seur@meta.data$celltype1 %>% table %>% data.frame %>% set_colnames(c("CellTypes","Number")) %>% dplyr::mutate(Per=100*Number/sum(Number))
})
names(EachSampleClusterDisB) <- unique(WholeTissueData_B_ident@meta.data$orig.ident)
EachSampleClusterDisC <- dplyr::bind_rows(EachSampleClusterDisB) %>% dplyr::mutate(Samples=rep(names(EachSampleClusterDisB),times=unlist(lapply(EachSampleClusterDisB,nrow))))
EachSampleClusterDisC$Tissue <- gsub("N.","NI",EachSampleClusterDisC$Samples)
EachSampleClusterDisC$Tissue <- gsub("P.*","AIP",EachSampleClusterDisC$Tissue)
pdf('Figure2c&SF2e.pdf',width = 6.8,height = 4)
ggplot(data = EachSampleClusterDisC,aes(x=Tissue,color=Tissue,y=Per))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1,size=1.5)+
  facet_wrap(~CellTypes,nrow =2,scales = 'free')+
  scale_x_discrete(limits = c('NI',"AIP"))+
  scale_y_continuous(expand = c(0.2,0))+
  scale_color_manual(limits=c("AIP","NI"),values=c("#CB181D","#1f78b4","#fdb462"))+
  ylab('Percentage(%)')+theme(axis.text = element_text(color = 'black'),panel.background = element_blank(),panel.grid = element_blank(),
                              legend.title = element_blank(),axis.text.y = element_text(color = 'black'),strip.background = element_blank(),
                              strip.text = element_text(size = 12),axis.title.x = element_blank(),axis.line = element_line(colour = 'black'),
                              legend.position = "none")+
  ggpubr::stat_compare_means(comparisons = list(c('NI','AIP')),paired = F,method = 'wilcox.test')
dev.off()

ClusterFreq <- WholeTissueData_B_ident@meta.data[,c("celltype1","Status")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("Naive B", "Memory B","ABC(IgD+)","ABC(IgD-)","GC B","Plasma G","Plasma A"))
pdf("SF2c.pdf",width = 10,height = 5)
ggplot(data = ClusterPer, aes(x = stim, fill=seurat_clusters,y=Per)) + 
  geom_bar(stat = "identity",width = 0.3)+
  scale_y_continuous(expand=c(0,0))+
  labs(y='percentage(%)')+
  coord_flip()+
  scale_fill_manual(values=.cluster_cols)+ 
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 17),axis.text.x=element_text(color = "black",size = 17),
        axis.title.y =element_blank(), axis.title.x =element_text(size = 20),axis.line = element_line(color="black"),legend.text = element_text(size=13))
dev.off()

#samples
ClusterFreq <- WholeTissueData_B_ident@meta.data[,c("celltype1","orig.ident")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("Naive B", "Memory B","ABC(IgD+)","ABC(IgD-)","GC B","Plasma G","Plasma A"))
ClusterPer$stim=factor(ClusterPer$stim,levels = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","N1","N2","N3","N4","N5","N6","N7","N8")) 
pdf("SF2d.pdf",width = 8,height = 6)
ggplot(data = ClusterPer, aes(x = stim, fill=seurat_clusters,y=Per)) + 
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=.cluster_cols)+ 
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 15),axis.text.x=element_text(color = "black",size = 15),
        axis.title.y = element_blank(),axis.title.x  = element_text(size = 18),axis.line = element_line(color="black"),legend.text = element_text(size=13))
dev.off()





