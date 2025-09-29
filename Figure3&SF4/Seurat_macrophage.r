setwd("/work/wjx/scRNA-seq/AIP/figure/figure3")
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
.cluster_cols <- c("#2E8B57","#FFA500","#6495ED","#DC050C","#00BFFF","#9400D3","#8B4513","#9400D3","#32CD32","#6A5ACD","#BDB76B","#8B4513","#DA70D6","#F8766D","#808080","#D3D3D3","#DAA520","#FFB6C1")


WholeTissueData_bigcluster=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure1/WholeTissueData_ident.rds.gz")
WholeTissueData_myeloid <- subset(WholeTissueData_bigcluster,celltype %in% c("Macrophage"))



WholeTissueData_myeloid  <- NormalizeData(WholeTissueData_myeloid) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_myeloid)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_myeloid <- ScaleData(WholeTissueData_myeloid,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_myeloid,ndims = 50)

WholeTissueData_myeloid <- FindNeighbors(WholeTissueData_myeloid,dims = 1:20)
WholeTissueData_myeloid  <- RunUMAP(object = WholeTissueData_myeloid, dims = 1:20) 
WholeTissueData_myeloid  <- FindClusters(WholeTissueData_myeloid, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_myeloid@meta.data,prefix = "RNA_snn_res.")

Idents(WholeTissueData_myeloid) <- "RNA_snn_res.0.4"



#remove doublets
WholeTissueData_myeloid_new <- subset(WholeTissueData_myeloid,RNA_snn_res.0.4!=4&RNA_snn_res.0.4!=5&RNA_snn_res.0.4!=6&RNA_snn_res.0.4!=7&RNA_snn_res.0.4!=8)
WholeTissueData_myeloid_new  <- NormalizeData(WholeTissueData_myeloid_new) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_myeloid_new)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_myeloid_new <- ScaleData(WholeTissueData_myeloid_new,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_myeloid_new,ndims = 50)

WholeTissueData_myeloid_new <- FindNeighbors(WholeTissueData_myeloid_new,dims = 1:15)
WholeTissueData_myeloid_new  <- RunUMAP(object = WholeTissueData_myeloid_new, dims = 1:15) 
WholeTissueData_myeloid_new  <- FindClusters(WholeTissueData_myeloid_new, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)


clustree(WholeTissueData_myeloid_new@meta.data,prefix = "RNA_snn_res.")


WholeTissueData_myeloid_new1 <- subset(WholeTissueData_myeloid_new,RNA_snn_res.0.9!=7&RNA_snn_res.0.9!=1&RNA_snn_res.0.9!=3)

WholeTissueData_myeloid_new1  <- NormalizeData(WholeTissueData_myeloid_new1) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_myeloid_new1)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_myeloid_new1 <- ScaleData(WholeTissueData_myeloid_new1,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_myeloid_new1,ndims = 50)

WholeTissueData_myeloid_new1 <- FindNeighbors(WholeTissueData_myeloid_new1,dims = 1:15)
WholeTissueData_myeloid_new1  <- RunUMAP(object = WholeTissueData_myeloid_new1, dims = 1:15) 
WholeTissueData_myeloid_new1  <- FindClusters(WholeTissueData_myeloid_new1, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)


clustree(WholeTissueData_myeloid_new1@meta.data,prefix = "RNA_snn_res.")

table(WholeTissueData_myeloid_new1@meta.data$RNA_snn_res.0.8)
table(WholeTissueData_myeloid_new1@meta.data$RNA_snn_res.0.7,WholeTissueData_myeloid_new1@meta.data$orig.ident)



levels(Idents(WholeTissueData_myeloid_new1))
Idents(WholeTissueData_myeloid_new1) <- 'RNA_snn_res.0.7'
WholeTissueData_myeloid_new1_res0.7.markers <- FindAllMarkers(WholeTissueData_myeloid_new1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)


.cluster_cols <- c("#DC050C","#FF7F00","#33A02C","#1965B0","#B17BA6")

WholeTissueData_myeloid_ident <- RenameIdents(WholeTissueData_myeloid_new1, "3"= "CXCL9+ Macrophage","2"= "CXCL13+ Macrophage","1"="TREM2+ Macrophage", "0"="FOLR2+ Macrophage" ,"5"= "FOLR2+ Macrophage","4"= "FCN1+ Macrophage")


WholeTissueData_myeloid_ident$celltype1=WholeTissueData_myeloid_ident@active.ident
table(WholeTissueData_myeloid_ident$celltype1,WholeTissueData_myeloid_ident$orig.ident)
table(WholeTissueData_myeloid_ident$Status)


readr::write_rds(WholeTissueData_myeloid_ident,path="WholeTissueData_macrophage_ident.rds.gz",compress="gz")



WholeTissueData_myeloid_ident <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure5/WholeTissueData_macrophage_ident.rds.gz")

.cluster_cols <- c("#DC050C","#FF7F00","#33A02C","#1965B0","#B17BA6")

WholeTissueData_myeloid_ident$Status <- factor(WholeTissueData_myeloid_ident$Status,levels = c("Normal","AIP"))

tiff('Figure3a.tiff',width = 1400,height = 700,res=150)
DimPlot(WholeTissueData_myeloid_ident,reduction = 'umap', label = F,repel = T,cols = .cluster_cols,group.by = "celltype1",split.by = "Status",pt.size = 0.5)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

WholeTissueData_myeloid_ident$orig.ident <- factor(WholeTissueData_myeloid_ident$orig.ident,levels = c("N1","N2","N3","N4","N5","N6","N7","N8","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))
tiff('SF4b.tiff',width = 1000,height = 910,res = 150)
DimPlot(WholeTissueData_myeloid_ident,reduction = 'umap', label = F,repel = T, group.by = "orig.ident")+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


tiff("SF4c.tiff",width = 1650,height = 1000,res =200)
FeaturePlot(WholeTissueData_myeloid_ident,c("CD68","CXCL9","CXCL13","TREM2","FOLR2","FCN1"),cols = rev(colorRampPalette(c(RColorBrewer::brewer.pal(10,'Spectral')[2:10]),space='rgb')(20)),ncol = 3)&
  theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line =element_blank())
dev.off()


Idents(WholeTissueData_myeloid_ident) <- "celltype1"
myeloid_markers <- FindAllMarkers(WholeTissueData_myeloid_ident,only.pos = T,min.pct = 0.1,logfc.threshold = 0.25)
write.csv(myeloid_markers,"/work/wjx/scRNA-seq/AIP/figure/figure5/myeloid_markers.csv")



PlotGenes <- c("MMP9","UBD","PLA2G2D","PTGDS","NR1H3","SPARC","CHIT1","IL32","ATOX1","CXCL9",
               "LBH","PPM1N","GNG2","CD5L","TMSB4X","FEZ1","HCAR3","UCP2","EBI3", "CXCL13",
               "OLFML3","HTRA1","GLDN","HPGDS","CCND1","ALOX5AP","PLD4","SPP1","TREM2","AXL",
               "CNRIP1", "SLC40A1","MS4A4A", "CTSF","SESN1","SELENOP","HLA-DPB1","CD163","FOLR2","PLD3",
               "VCAN","CD52","S100A9","FCN1","S100A8","FLNA","EMP3","S100A10","S100A6","S100A4")

Sub_Data <- WholeTissueData_myeloid_ident

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

pdf("Figure3b.pdf",width = 8,height = 10)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),
                        show_row_names = F,show_column_names = T  ,column_names_rot = 45,top_annotation = Sub_ha_bar)
dev.off()





EachSampleClusterDisB <- lapply(unique(WholeTissueData_myeloid_ident@meta.data$orig.ident),function(i){
  seur <- subset(WholeTissueData_myeloid_ident,cells=rownames(WholeTissueData_myeloid_ident@meta.data[WholeTissueData_myeloid_ident@meta.data$orig.ident %in% i,]))
  seur@meta.data$celltype1 %>% table %>% data.frame %>% set_colnames(c("CellTypes","Number")) %>% dplyr::mutate(Per=100*Number/sum(Number))
})
names(EachSampleClusterDisB) <- unique(WholeTissueData_myeloid_ident@meta.data$orig.ident)
EachSampleClusterDisC <- dplyr::bind_rows(EachSampleClusterDisB) %>% dplyr::mutate(Samples=rep(names(EachSampleClusterDisB),times=unlist(lapply(EachSampleClusterDisB,nrow))))
EachSampleClusterDisC$Tissue <- gsub("N.","NI",EachSampleClusterDisC$Samples)
EachSampleClusterDisC$Tissue <- gsub("P.*","AIP",EachSampleClusterDisC$Tissue)
pdf('Figure3c&SF4f.pdf',width = 6.8,height = 4)
ggplot(data = EachSampleClusterDisC,aes(x=Tissue,color=Tissue,y=Per))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1,size=1.5)+
  facet_wrap(~CellTypes,ncol =4,scales = 'free')+
  scale_x_discrete(limits = c('NI',"AIP"))+
  scale_y_continuous(expand = c(0.2,0))+
  scale_color_manual(limits=c("AIP","NI"),values=c("#CB181D","#1f78b4","#fdb462"))+
  ylab('Percentage(%)')+theme(axis.text = element_text(color = 'black'),panel.background = element_blank(),panel.grid = element_blank(),
                              legend.title = element_blank(),axis.text.y = element_text(color = 'black'),strip.background = element_blank(),
                              strip.text = element_text(size = 12),axis.title.x = element_blank(),axis.line = element_line(colour = 'black'),
                              legend.position = "none")+
  ggpubr::stat_compare_means(comparisons = list(c('NI','AIP')),paired = F,method = 'wilcox.test')
dev.off()

ClusterFreq <- WholeTissueData_myeloid_ident@meta.data[,c("celltype1","Status")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("CXCL9+ Macrophage","CXCL13+ Macrophage","TREM2+ Macrophage","FOLR2+ Macrophage","FCN1+ Macrophage"))
ClusterPer$stim=factor(ClusterPer$stim,levels = c("AIP","Normal")) 
pdf("SF4d.pdf",width = 10,height = 5)
ggplot(data = ClusterPer, aes(x = stim, fill=seurat_clusters,y=Per)) + #log2(as.numeric(CopyNumber)))
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


ClusterFreq <- WholeTissueData_myeloid_ident@meta.data[,c("celltype1","orig.ident")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("CXCL9+ Macrophage","CXCL13+ Macrophage","TREM2+ Macrophage","FOLR2+ Macrophage","FCN1+ Macrophage"))
ClusterPer$stim=factor(ClusterPer$stim,levels = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","N3","N4","N5","N6","N7","N8")) 
ClusterPer <- ClusterPer[!is.na(ClusterPer$stim),]
pdf("SF4e.pdf",width = 8,height = 6)
ggplot(data = ClusterPer, aes(x = stim, fill=seurat_clusters,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=.cluster_cols)+ 
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 15),axis.text.x=element_text(color = "black",size = 15),
        axis.title.y = element_blank(),axis.title.x  = element_text(size = 18),axis.line = element_line(color="black"),legend.text = element_text(size=13))
dev.off()
