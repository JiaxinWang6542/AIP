setwd("/work/wjx/scRNA-seq/AIP")
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
options(stringsAsFactors = F)
set.seed(42)

.cluster_cols <- c(
   "#FB8072","#DC050C", "#7BAFDE","#1965B0", "#B2DF8A","#33A02C",  "#8DD3C7","#E7298A", "#E78AC3",
  "#882E72","#B17BA6", "#FF7F00", "#FDB462", "#A6761D","#E6AB02", "#7570B3", "#BEAED4", "#666666", 
  "#999999","#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000","#1e90ff", "#00bfff", "#56ff0d","#aeae5c",  "#ffff00")



##All_Tissue_samples
PatientList <- grep(pattern = "Inflam|Normal",list.dirs(),value=T)
PatientList
WholeTissue_object_list <- lapply(PatientList,function(x){
  Xdata <- Read10X(data.dir = paste( x,"/",sep=""))
  XF <- CreateSeuratObject(counts =  Xdata, project = x,min.cells = 3)
  double_score=read.csv(paste(paste("doublets",gsub("./matrix/Tissue","/doublets_result",x),sep=""),"_scrublet_doubletdetection_mark_doublets_barcodes.csv",sep=""))
  XF$scrublet_doublet_scores=double_score$scrublet_doublet_scores
  XF$scrublet_doublets=double_score$scrublet_doublets
  XF$doubletdetection_doublets=double_score$doubletdetection_doublets
  print(x)
  return(XF)
})

PatientList
names(WholeTissue_object_list) <- gsub("./matrix/Tissue/","",PatientList)
names(WholeTissue_object_list)
WholeTissueData <- WholeTissue_object_list[[1]]
WholeTissueData <- merge(x = WholeTissueData, y = WholeTissue_object_list[2:length(WholeTissue_object_list)],add.cell.ids=c(names(WholeTissue_object_list)),merge.data=T,project = "SeuratProject")
head(rownames(WholeTissueData@meta.data))
WholeTissueData@meta.data$Sample <-   gsub("_Inflam.*|_Normal.*","",rownames(WholeTissueData@meta.data))
table(WholeTissueData@meta.data$Sample)

head(WholeTissueData@meta.data)
WholeTissueData$orig.ident=gsub("./matrix/Tissue/","",WholeTissueData$orig.ident)
WholeTissueData$orig.ident=gsub("_AIP_.*","",WholeTissueData$orig.ident)
WholeTissueData$orig.ident=gsub("P1_1|P1_2","P1",WholeTissueData$orig.ident)
table(WholeTissueData@meta.data$orig.ident)
head(WholeTissueData@meta.data)
mito.genes <- grep(pattern = "^MT-", x = rownames(WholeTissueData@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(WholeTissueData@assays[["RNA"]][mito.genes, ])/Matrix::colSums(WholeTissueData@assays[["RNA"]])

WholeTissueData <- AddMetaData(object = WholeTissueData, metadata = percent.mito, col.name = "percent.mito") 

head(WholeTissueData@meta.data)
readr::write_rds(WholeTissueData,path="/work/wjx/scRNA-seq/AIP/figure/figure1/WholeTissueData_Raw.rds.gz",compress="gz")
rm(WholeTissue_object_list,mito.genes,PatientList,percent.mito)

WholeTissueData=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure1/WholeTissueData_Raw.rds.gz")
setwd("/work/wjx/scRNA-seq/AIP/figure/figure1")

VlnPlot(object =WholeTissueData, features = c("percent.mito"),group.by="orig.ident",pt.size = 0)+
  geom_hline(yintercept = seq(0.05,0.5,by=0.05),linetype="dashed")


VlnPlot(object = WholeTissueData, features = c("nFeature_RNA"),group.by="orig.ident",pt.size = 0)+
  geom_hline(yintercept = seq(1000,8000,by=1000),linetype="dashed",color="red")


VlnPlot(object = WholeTissueData, features = c("nCount_RNA"),group.by="orig.ident",pt.size = 0,y.max = 80000)+
  geom_hline(yintercept = seq(10000,100000,by=10000),linetype="dashed",color="red")


FeatureScatter(object = WholeTissueData, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident")



###Ribo genes
Ribo.genes <- grep(pattern = "^RPL|^RPS", x = rownames(WholeTissueData@assays[["RNA"]]), value = TRUE)
WholeTissueData$percent.ribo <- Matrix::colSums(WholeTissueData@assays[["RNA"]][Ribo.genes, ])/Matrix::colSums(WholeTissueData@assays[["RNA"]])
VlnPlot(object =WholeTissueData, features = c("percent.ribo"),group.by="orig.ident",pt.size = 0)



###blood genes
HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
HB_m <- match(HB.genes_total,rownames(WholeTissueData@assays$RNA))
HB.genes <- rownames(WholeTissueData@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
WholeTissueData[["percent.HB"]]<-PercentageFeatureSet(WholeTissueData,features=HB.genes)
VlnPlot(object =WholeTissueData, features = c("percent.HB"),group.by="orig.ident")+
  geom_hline(yintercept = seq(10,100,by=10),linetype="dashed")





###cutoff
WholeTissueData_cutoff=subset(WholeTissueData,subset=scrublet_doublets == "False")###remove doublets
WholeTissueData_cutoff=subset(WholeTissueData_cutoff,subset=doubletdetection_doublets == "False")###remove doublets
WholeTissueData_cutoff <- subset(WholeTissueData_cutoff, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mito >  -Inf & percent.mito < 0.1& percent.HB < 25 )


##pca-umap clustering
WholeTissueData_log  <- NormalizeData(WholeTissueData_cutoff) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_log)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_log <- ScaleData(WholeTissueData_log,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissueData_log,ndims = 50)
WholeTissueData_log <- FindNeighbors(WholeTissueData_log,dims = 1:20)
WholeTissueData_log  <- RunUMAP(object = WholeTissueData_log, dims = 1:20) 
WholeTissueData_log   <- FindClusters(WholeTissueData_log , resolution = c(seq(0.1,0.4,0.1)), algorithm = 1)

clustree(WholeTissueData_log@meta.data,prefix = "RNA_snn_res.")

levels(Idents(WholeTissueData_log))
Idents(WholeTissueData_log) <- 'RNA_snn_res.0.2'
library(future)
plan(multisession,workers=20)
options(future.globals.maxSize= 2097152000)
WholeTissueData_log_markers <- FindMarkers(WholeTissueData_log, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


##remove 7(platelet)、stamoch epithelial(15 17)、enterocyte (16 18)
WholeTissueData_new <- subset(WholeTissueData_log, subset = RNA_snn_res.0.2!=7 &  RNA_snn_res.0.2!=15 & RNA_snn_res.0.2!=16&  RNA_snn_res.0.2!=17 & RNA_snn_res.0.2!=18)
plan(sequential)

##pca-umap clustering
#wihtout harmony
WholeTissueData_new  <- NormalizeData(WholeTissueData_new) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_new)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_new <- ScaleData(WholeTissueData_new,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_new,ndims = 50)

WholeTissueData_new <- FindNeighbors(WholeTissueData_new,dims = 1:20)
WholeTissueData_new  <- RunUMAP(object = WholeTissueData_new, dims = 1:20) 
WholeTissueData_new   <- FindClusters(WholeTissueData_new , resolution = c(seq(0.1,0.4,0.1)), algorithm = 1)

clustree(WholeTissueData_new@meta.data,prefix = "RNA_snn_res.")

Idents(WholeTissueData_new) <- 'RNA_snn_res.0.3'
library(future)
plan(multisession,workers=20)
options(future.globals.maxSize= 2097152000)


WholeTissueData_res0.3_new.markers <- FindAllMarkers(WholeTissueData_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



#Identify celltypes
WholeTissueData_ident <- RenameIdents(WholeTissueData_new, "0" = "T","2" = "T","3" = "T","11" = "T","14" = "T","18" = "T","5" = "NK","4"= "B","7"="Plasma","1"="Monocyte","8"="Monocyte","9"="Macrophage","17"="DC","10"="Neutrophil","16"="Prolifer","6"="Ductal","19"="Ductal","13"="Acinar","12"="Endothelial","15"="Fibroblast","20"="Endocrine")
WholeTissueData_ident$celltype=WholeTissueData_ident@active.ident
WholeTissueData_ident$Status=gsub("N.*","Normal",gsub("P.*","AIP",WholeTissueData_ident$orig.ident))
WholeTissueData_ident_markers <- FindAllMarkers(WholeTissueData_ident,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(WholeTissueData_ident_markers ,file = "/work/wjx/scRNA-seq/AIP/figure/figure1/bigcluster_markers_final.csv")
WholeTissueData_ident_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_logFC) -> WholeTissueData_ident_markers_top20
readr::write_rds(WholeTissueData_ident,path="./result/WholeTissueData_ident.rds.gz",compress="gz")
WholeTissueData_ident <- readRDS("./result/WholeTissueData_ident.rds.gz")

##plot
table(WholeTissueData_ident$Status)
plot1 <- DimPlot(subset(WholeTissueData_ident,Status=="AIP"),pt.size=0.2,label = F,cols = .cluster_cols)+
  labs(title="AIP(47586 cells)")+
  theme(legend.text=element_text(size=15),legend.key.size = unit(20,"pt"),plot.title = element_text(hjust = 0.5,size=20))
plot2 <- DimPlot(subset(WholeTissueData_ident,Status=="Normal"),pt.size=0.2,label = F,cols = .cluster_cols)+
  labs(title="Normal(31793 cells)")+
  theme(legend.text=element_text(size=15),legend.key.size = unit(20,"pt"),plot.title = element_text(hjust = 0.5,size=20)) 

tiff('Figure1b.tiff',width = 3500,height = 1400,res = 250)
plot2+plot1+ plot_layout(guides = 'collect')
dev.off()

WholeTissueData_ident$orig.ident <- factor(WholeTissueData_ident$orig.ident,levels = c("N1","N2","N3","N4","N5","N6","N7","N8","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))
tiff('SF1a.tiff',width = 1120,height = 910,res = 150)
DimPlot(WholeTissueData_ident,reduction = 'umap', label = F,repel = T, group.by = "orig.ident")+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


tiff("SF1b.tiff",width = 2500,height = 1500,res =200)
FeaturePlot(WholeTissueData_ident,c("CD4","CD8A","NCAM1","MS4A1","MZB1","CD14","C1QC","FCER1A","FCGR3B","MKI67","EPCAM","PRSS1","PECAM1","COL1A1","TTR"),cols = rev(colorRampPalette(c(RColorBrewer::brewer.pal(10,'Spectral')[2:10]),space='rgb')(20)),ncol = 5)&
  theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line =element_blank())
dev.off()



ClusterFreq <- WholeTissueData_ident@meta.data[,c("celltype","Status")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("T","NK","B","Plasma","Monocyte","Macrophage","DC","Neutrophil","Prolifer","Ductal","Acinar","Endothelial","Fibroblast","Endocrine"))
ClusterPer$stim=factor(ClusterPer$stim,levels = c("AIP","Normal")) 
pdf("Figure1d.pdf",width = 10,height = 5)
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


EachSampleClusterDisB <- lapply(unique(WholeTissueData_ident@meta.data$orig.ident),function(i){
  seur <- subset(WholeTissueData_ident,cells=rownames(WholeTissueData_ident@meta.data[WholeTissueData_ident@meta.data$orig.ident %in% i,]))
  seur@meta.data$celltype %>% table %>% data.frame %>% set_colnames(c("CellTypes","Number")) %>% dplyr::mutate(Per=100*Number/sum(Number))
})
names(EachSampleClusterDisB) <- unique(WholeTissueData_ident@meta.data$orig.ident)
EachSampleClusterDisC <- dplyr::bind_rows(EachSampleClusterDisB) %>% dplyr::mutate(Samples=rep(names(EachSampleClusterDisB),times=unlist(lapply(EachSampleClusterDisB,nrow))))
EachSampleClusterDisC$Tissue <- gsub("N.","NI",EachSampleClusterDisC$Samples)
EachSampleClusterDisC$Tissue <- gsub("P.*","AIP",EachSampleClusterDisC$Tissue)
pdf('Figure1e&SF1d.pdf.pdf',width = 8.5,height = 6)
ggplot(data = EachSampleClusterDisC,aes(x=Tissue,color=Tissue,y=Per))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1,size=1.5)+
  facet_wrap(~CellTypes,nrow =3,scales = 'free')+
  scale_x_discrete(limits = c('NI',"AIP"))+
  scale_y_continuous(expand = c(0.2,0))+
  scale_color_manual(limits=c("AIP","NI"),values=c("#CB181D","#1f78b4","#fdb462"))+
  ylab('Percentage(%)')+theme(axis.text = element_text(color = 'black'),panel.background = element_blank(),panel.grid = element_blank(),
                              legend.title = element_blank(),axis.text.y = element_text(color = 'black'),strip.background = element_blank(),
                              strip.text = element_text(size = 12),axis.title.x = element_blank(),axis.line = element_line(colour = 'black'),
                              legend.position = "none")+
  ggpubr::stat_compare_means(comparisons = list(c('NI','AIP')),paired = F,method = 'wilcox.test')
dev.off()


#samples
ClusterFreq <- WholeTissueData_ident@meta.data[,c("celltype","orig.ident")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("T","NK","B","Plasma","Monocyte","Macrophage","DC","Neutrophil","Prolifer","Ductal","Acinar","Endothelial","Fibroblast","Endocrine"))
ClusterPer$stim=factor(ClusterPer$stim,levels = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","N1","N2","N3","N4","N5","N6","N7","N8")) 
pdf("SF1c.pdf",width = 8,height = 6)
ggplot(data = ClusterPer, aes(x = stim, fill=seurat_clusters,y=Per)) + 
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=.cluster_cols)+ 
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 15),axis.text.x=element_text(color = "black",size = 15),
        axis.title.y = element_blank(),axis.title.x  = element_text(size = 18),axis.line = element_line(color="black"),legend.text = element_text(size=13))
dev.off()



PlotGenes <- WholeTissueData_ident_markers_top20$gene
Sub_Data <- WholeTissueData_ident
Sub_TF_Matrix <- Sub_Data@assays$RNA@data[PlotGenes,] %>% as.matrix()
Sub_TF_Matrix <- Sub_TF_Matrix[apply(Sub_TF_Matrix,1,sum) > 0,]
Sub_ClusterID <- lapply(split(Sub_Data@meta.data,list(Sub_Data$celltype)), function(x)rownames(x))
Sub_TF_ClusterMean <- t(apply(Sub_TF_Matrix,1,function(x){
  lapply(Sub_ClusterID,function(ID){
    mean(x[ID])
  }) %>% unlist
}))
Sub_TF_ClusterMean <- Sub_TF_ClusterMean[apply(Sub_TF_ClusterMean,1,sum)>0,]
Sub_TF_ClusterMean[,1:ncol(Sub_TF_ClusterMean)] <- t(apply(Sub_TF_ClusterMean,1,scale))

Sub_barcolor <- .cluster_cols[1:14]
names(Sub_barcolor) <- unique(Sub_Data$celltype)
Sub_SampleBar <-  unique(Sub_Data$celltype)
names(Sub_SampleBar) <-  unique(Sub_Data$celltype)

library(ComplexHeatmap)
Sub_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Sub_SampleBar,
                                                col = list(bar = Sub_barcolor),show_legend = F,show_annotation_name = F)

pdf("Figure1c.pdf",width = 10.5,height = 8)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),
                        show_row_names = F,show_column_names = T  ,column_names_rot = 45,top_annotation = Sub_ha_bar)
dev.off()
