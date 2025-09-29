setwd("/work/wjx/scRNA-seq/AIP/figure/figure4")
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

######T cell
table(WholeTissueData_ident$celltype1)
WholeTissueData_T <- subset(WholeTissueData_ident,celltype1 %in% c("T"))
table(WholeTissueData_T@meta.data$celltype1)
WholeTissueData_T  <- NormalizeData(WholeTissueData_T) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_T)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_T <- ScaleData(WholeTissueData_T,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissueData_T,ndims = 50)
WholeTissueData_T <- FindNeighbors(WholeTissueData_T,dims = 1:20)
WholeTissueData_T  <- RunUMAP(object = WholeTissueData_T, dims = 1:20) 
WholeTissueData_T  <- FindClusters(WholeTissueData_T, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_T@meta.data,prefix = "RNA_snn_res.")



Idents(WholeTissueData_T) <- 'RNA_snn_res.0.9'
library(future)
plan(multisession,workers=30)
options(future.globals.maxSize= 2097152000)
WholeTissueData_T_res0.9.markers <- FindAllMarkers(WholeTissueData_T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



###remove doublets
WholeTissueData_T_new <- subset(WholeTissueData_T,RNA_snn_res.0.9!=8 & RNA_snn_res.0.9!=11& RNA_snn_res.0.9!=13& RNA_snn_res.0.9!=16& RNA_snn_res.0.9!=17& RNA_snn_res.0.9!=18& RNA_snn_res.0.9!=20)
plan(sequential)
WholeTissueData_T_new  <- NormalizeData(WholeTissueData_T_new) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_T_new)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_T_new <- ScaleData(WholeTissueData_T_new,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissueData_T_new,ndims = 50)
WholeTissueData_T_new <- FindNeighbors(WholeTissueData_T_new,dims = 1:20)
WholeTissueData_T_new  <- RunUMAP(object = WholeTissueData_T_new, dims = 1:20) 
WholeTissueData_T_new  <- FindClusters(WholeTissueData_T_new, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_T_new@meta.data,prefix = "RNA_snn_res.")
Idents(WholeTissueData_T_new) <- "RNA_snn_res.1.2"


###remove doublets
WholeTissueData_T_new1 <- subset(WholeTissueData_T_new,RNA_snn_res.1.2!=13&RNA_snn_res.1.2!=14&RNA_snn_res.1.2!=17&RNA_snn_res.1.2!=18)
plan(sequential)
WholeTissueData_T_new1  <- NormalizeData(WholeTissueData_T_new1) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_T_new1)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_T_new1 <- ScaleData(WholeTissueData_T_new1,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_T_new1,ndims = 50)

WholeTissueData_T_new1 <- FindNeighbors(WholeTissueData_T_new1,dims = 1:20)
WholeTissueData_T_new1  <- RunUMAP(object = WholeTissueData_T_new1, dims = 1:20) 
WholeTissueData_T_new1  <- FindClusters(WholeTissueData_T_new1, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_T_new1@meta.data,prefix = "RNA_snn_res.")



Idents(WholeTissueData_T_new1) <- 'RNA_snn_res.1.4'
library(future)
plan(multisession,workers=20)
options(future.globals.maxSize= 2097152000)
WholeTissueData_T_new1_res1.4markers <- FindAllMarkers(WholeTissueData_T_new1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


WholeTissueData_T_ident <- RenameIdents(WholeTissueData_T_new1, "3" = "CD4+ Tn","6" = "CD4+ Tn","9" = "CD4+ Tn","17" = "CD4+ Tn","4"="CD8+ Tn","0" = "LIMS1+ Tm","1" = "TIMP1+ Tm","5" = "TIMP1+ Tm","13" = "TIMP1+ Tm","2" = "Tem","11" = "Tem","10"= "CD4+ Temra","8"= "CD8+ Temra","15"= "Trm","12" = "Tfh","7"= "Treg","18"= "Treg","14"= "MAIT","16"= "gdT","19"= "ISG+ T")

WholeTissueData_T_ident$celltype1=WholeTissueData_T_ident@active.ident


.cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A","#A6761D", "#55A1B1", "#8DD3C7", 
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

WholeTissueData_T_ident$Status <- factor(WholeTissueData_T_ident$Status,levels = c("Normal","AIP"))

tiff('SF5a.tiff',width = 1400,height = 700,res=150)
DimPlot(WholeTissueData_T_ident,reduction = 'umap', label = F,repel = T,cols = .cluster_cols,group.by = "celltype1",split.by = "Status",pt.size = 0.5)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf("SF5b.pdf",width = 10,height = 4.5)
Seurat::DotPlot(WholeTissueData_T_ident,
                features = c("CD3D","CD4","CD8A","CCR7","SELL","IL7R","S100A4","LIMS1","TIMP1","GZMK","CCL5","GZMH","FGFBP2","ITGAE","ZNF683","PDCD1","CXCL13","FOXP3","IL2RA","TNFRSF18","SLC4A10","TRDV2","TRGV9","IFIT1","ISG15"),
                cols = 'RdYlBu')+
  theme(axis.text.x  = element_text(angle = 45,vjust = .9,hjust = 0.8))
dev.off()

readr::write_rds(WholeTissueData_T_ident,path="result/WholeTissueData_T_ident",compress="gz")



###CD4 T
WholeTissueData_CD4T_raw <- subset(WholeTissueData_T_ident,celltype1 %in% c("CD4+ Tn","LIMS1+ Tm","TIMP1+ Tm","CD4+ Temra","Tfh","Treg","ISG+ T"))
table(WholeTissueData_CD4T_raw@meta.data$celltype1)

plan(sequential)
WholeTissueData_CD4T_raw  <- NormalizeData(WholeTissueData_CD4T_raw) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_CD4T_raw)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_CD4T_raw <- ScaleData(WholeTissueData_CD4T_raw,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_CD4T_raw,ndims = 50)

WholeTissueData_CD4T_raw <- FindNeighbors(WholeTissueData_CD4T_raw,dims = 1:20)
WholeTissueData_CD4T_raw  <- RunUMAP(object = WholeTissueData_CD4T_raw, dims = 1:20) 
WholeTissueData_CD4T_raw  <- FindClusters(WholeTissueData_CD4T_raw, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_CD4T_raw@meta.data,prefix = "RNA_snn_res.")


Idents(WholeTissueData_CD4T_raw) <- 'RNA_snn_res.0.9'


#remove doublets
WholeTissueData_CD4T <- subset(WholeTissueData_CD4T_raw,RNA_snn_res.0.9!=10 & RNA_snn_res.0.9!=16 )

WholeTissueData_CD4T  <- NormalizeData(WholeTissueData_CD4T) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_CD4T)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_CD4T <- ScaleData(WholeTissueData_CD4T,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissueData_CD4T,ndims = 50)
WholeTissueData_CD4T <- FindNeighbors(WholeTissueData_CD4T,dims = 1:15)
WholeTissueData_CD4T  <- RunUMAP(object = WholeTissueData_CD4T, dims = 1:15) 
WholeTissueData_CD4T  <- FindClusters(WholeTissueData_CD4T, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_CD4T@meta.data,prefix = "RNA_snn_res.")


levels(Idents(WholeTissueData_CD4T))
Idents(WholeTissueData_CD4T) <- 'RNA_snn_res.1.2'
WholeTissueData_CD4T_res1.2markers <- FindAllMarkers(WholeTissueData_CD4T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



WholeTissueData_CD4T_ident <- RenameIdents(WholeTissueData_CD4T, "2" = "Tn","5" = "Tn","9" = "Tn","10" = "Tn","11" = "Tn","0" = "LIMS1+ Tm","8" = "LIMS1+ Tm","3" = "TIMP1+ Tm","6" = "TIMP1+ Tm","14" = "TIMP1+ Tm","17" = "TIMP1+ Tm","1" = "Tem","7" = "Temra","4"= "TNFRSF18- Treg","13"= "TNFRSF18+ Treg","12" = "Tfh","15" = "Tfh","16" = "ISG+ T")



table(WholeTissueData_CD4T_ident$celltype1)
WholeTissueData_CD4T_ident$celltype2=WholeTissueData_CD4T_ident@active.ident



########CD4T
.cluster_cols <- c(
  "#FB8072", "#1965B0", "#7BAFDE",  "#33A02C", "#B2DF8A","#FF7F00", "#FDB462","#DC050C","#882E72")

WholeTissueData_CD4T_ident$Status <- factor(WholeTissueData_CD4T_ident$Status,levels = c("Normal","AIP"))

tiff('Figure4a.tiff',width = 1400,height = 700,res=150)
DimPlot(WholeTissueData_CD4T_ident,reduction = 'umap', label = F,repel = T,cols = .cluster_cols,group.by = "celltype2",split.by = "Status",pt.size = 0.5)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

WholeTissueData_CD4T_ident$orig.ident <- factor(WholeTissueData_CD4T_ident$orig.ident,levels = c("N1","N2","N3","N4","N5","N6","N7","N8","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"))

tiff('SF5c.tiff',width = 1120,height = 910,res = 150)
DimPlot(WholeTissueData_CD4T_ident,reduction = 'umap', label = F,repel = T, group.by = "orig.ident")+#NoLegend()+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

tiff("SF5d.tiff",width = 2200,height = 1500,res =200)
FeaturePlot(WholeTissueData_CD4T_ident,c("CD3D","CD4","CCR7","ANXA1","LIMS1","TIMP1","GZMK","GZMB","FOXP3","TNFRSF18","PDCD1","ISG15"),cols = rev(colorRampPalette(c(RColorBrewer::brewer.pal(10,'Spectral')[2:10]),space='rgb')(20)),ncol = 4)&
  theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line =element_blank(), legend.position = 'none')
dev.off()

Idents(WholeTissueData_CD4T_ident) <- "celltype2"
CD4T_markers <- FindAllMarkers(WholeTissueData_CD4T_ident,only.pos = T,min.pct = 0.1,logfc.threshold = 0.25)


PlotGenes <- c("TMIGD2","NOG","BACH2","CHRM3-AS2","AIF1","ACTN1","CCR7","LEF1","SELL","TCF7",
               "LIMS1","PASK","RNASET2","FYB1","RNASET2","CNN2","LINC00402","TRIB2","PRKCB","TPR",
               "LINC02762","FUT7","WDR86-AS1","CTSH","CAPG","LMNA","USP46","CCR6","PI16","TIMP1",
               "GZMK","IFNG-AS1","CCR2","MT-ND1","MT-ND2","NSG1","DPP4","CEBPD","CXCR3","IFNGR1",
               "NKG7","GNLY","GZMH","FGFBP2","GZMB","PRF1","CX3CR1","EFHD2","CCL4","GZMA",
               "CCR4","SHMT2","CCR10","HPGD","USP15","LFNG","RNF214","IL10RA","IL2RA","FOXP3",
               "CXCR6","LINC01943","BATF","GADD45A","TNFAIP3","CKLF","TNFRSF18","FOXP3","PTTG1","TNFRSF1B",
               "CXCL13","PDCD1","IL21","CD200","AC004585.1","GNG4","DRAIC","ITM2A","LINC01480","IGFBP4",
               "IFIT1","IFIT3","IFI44L","IFI6","ISG15","IFI44","IFI35","CMPK2","MX1","OAS3")

Sub_Data <- WholeTissueData_CD4T_ident

Sub_TF_Matrix <- Sub_Data@assays$RNA@data[PlotGenes,] %>% as.matrix()

Sub_TF_Matrix <- Sub_TF_Matrix[apply(Sub_TF_Matrix,1,sum) > 0,]

Sub_ClusterID <- lapply(split(Sub_Data@meta.data,list(Sub_Data$celltype2)), function(x)rownames(x))

Sub_TF_ClusterMean <- t(apply(Sub_TF_Matrix,1,function(x){
  lapply(Sub_ClusterID,function(ID){
    mean(x[ID])
  }) %>% unlist
}))

Sub_TF_ClusterMean <- Sub_TF_ClusterMean[apply(Sub_TF_ClusterMean,1,sum)>0,]
Sub_TF_ClusterMean[,1:ncol(Sub_TF_ClusterMean)] <- t(apply(Sub_TF_ClusterMean,1,scale))

Sub_barcolor <- .cluster_cols
names(Sub_barcolor) <- unique(Sub_Data$celltype2)
Sub_SampleBar <-  unique(Sub_Data$celltype2)
names(Sub_SampleBar) <-  unique(Sub_Data$celltype2)


library(ComplexHeatmap)
Sub_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Sub_SampleBar,
                                                col = list(bar = Sub_barcolor),show_legend = F,show_annotation_name = F)

pdf("Figure4b.pdf",width = 8,height = 10)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),
                        show_row_names = F,show_column_names = T  ,column_names_rot = 45,top_annotation = Sub_ha_bar)
dev.off()



EachSampleClusterDisB <- lapply(unique(WholeTissueData_CD4T_ident@meta.data$orig.ident),function(i){
  seur <- subset(WholeTissueData_CD4T_ident,cells=rownames(WholeTissueData_CD4T_ident@meta.data[WholeTissueData_CD4T_ident@meta.data$orig.ident %in% i,]))
  seur@meta.data$celltype2 %>% table %>% data.frame %>% set_colnames(c("CellTypes","Number")) %>% dplyr::mutate(Per=100*Number/sum(Number))
})
names(EachSampleClusterDisB) <- unique(WholeTissueData_CD4T_ident@meta.data$orig.ident)
EachSampleClusterDisC <- dplyr::bind_rows(EachSampleClusterDisB) %>% dplyr::mutate(Samples=rep(names(EachSampleClusterDisB),times=unlist(lapply(EachSampleClusterDisB,nrow))))
EachSampleClusterDisC$Tissue <- gsub("N.","NI",EachSampleClusterDisC$Samples)
EachSampleClusterDisC$Tissue <- gsub("P.*","AIP",EachSampleClusterDisC$Tissue)
pdf('Figure4c&SF5e.pdf',width = 6.8,height = 4)
ggplot(data = EachSampleClusterDisC,aes(x=Tissue,color=Tissue,y=Per))+#,group=Patients))+
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

ClusterFreq <- WholeTissueData_CD4T_ident@meta.data[,c("celltype2","Status")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("Tn", "LIMS1+ Tm","TIMP1+ Tm","Tem","Temra","TNFRSF18- Treg","TNFRSF18+ Treg","Tfh","ISG+ T"))
ClusterPer$stim=factor(ClusterPer$stim,levels = c("AIP","Normal")) 
pdf("SF5f.pdf",width = 10,height = 5)
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


ClusterFreq <- WholeTissueData_CD4T_ident@meta.data[,c("celltype2","orig.ident")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("Tn", "LIMS1+ Tm","TIMP1+ Tm","Tem","Temra","TNFRSF18- Treg","TNFRSF18+ Treg","Tfh","ISG+ T"))
ClusterPer$stim=factor(ClusterPer$stim,levels = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","N1","N2","N3","N4","N5","N6","N7","N8")) 
pdf("SF5g.pdf",width = 8,height = 6)
ggplot(data = ClusterPer, aes(x = stim, fill=seurat_clusters,y=Per)) + 
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=.cluster_cols)+ 
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black",size = 15),axis.text.x=element_text(color = "black",size = 15),
        axis.title.y = element_blank(),axis.title.x  = element_text(size = 18),axis.line = element_line(color="black"),legend.text = element_text(size=13))
dev.off()

readr::write_rds(WholeTissueData_CD4T_ident,path="result/WholeTissueData_CD4T_ident",compress="gz")



###CD8 
WholeTissueData_CD8T_raw <- subset(WholeTissueData_T_ident,celltype1 %in%c("CD8+ Tn","Tem","CD8+ Temra","Trm","MAIT") )
WholeTissueData_CD8T_raw  <- NormalizeData(WholeTissueData_CD8T_raw) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_CD8T_raw)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_CD8T_raw <- ScaleData(WholeTissueData_CD8T_raw,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_CD8T_raw,ndims = 50)

WholeTissueData_CD8T_raw <- FindNeighbors(WholeTissueData_CD8T_raw,dims = 1:15)
WholeTissueData_CD8T_raw  <- RunUMAP(object = WholeTissueData_CD8T_raw, dims = 1:15) 
WholeTissueData_CD8T_raw  <- FindClusters(WholeTissueData_CD8T_raw, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_CD8T_raw@meta.data,prefix = "RNA_snn_res.")


Idents(WholeTissueData_CD8T_raw) <- 'RNA_snn_res.0.8'


##remove doublets
WholeTissueData_CD8T <- subset(WholeTissueData_CD8T_raw,RNA_snn_res.0.8!=9 )
WholeTissueData_CD8T  <- NormalizeData(WholeTissueData_CD8T) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_CD8T)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_CD8T <- ScaleData(WholeTissueData_CD8T,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_CD8T,ndims = 50)

WholeTissueData_CD8T <- FindNeighbors(WholeTissueData_CD8T,dims = 1:15)
WholeTissueData_CD8T  <- RunUMAP(object = WholeTissueData_CD8T, dims = 1:15) 
WholeTissueData_CD8T  <- FindClusters(WholeTissueData_CD8T, resolution = c(seq(0.1,1.5,0.1)), algorithm = 1)

clustree(WholeTissueData_CD8T@meta.data,prefix = "RNA_snn_res.")

Idents(WholeTissueData_CD8T) <- 'RNA_snn_res.0.7'
WholeTissueData_CD8T_res0.7.markers <- FindAllMarkers(WholeTissueData_CD8T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



WholeTissueData_CD8T_ident <- RenameIdents(WholeTissueData_CD8T, "3" = "Tn","5" = "Tn","8" = "Tn","1" = "Tem","4" = "Tem","6" = "Tem","7"="Trm","11"="Trm","0"="Temra","9"="Temra","10"="Temra","12"="Temra","2"="MAIT")

WholeTissueData_CD8T_ident$celltype2=WholeTissueData_CD8T_ident@active.ident

.cluster_cols <- c(
  "#FB8072", "#1965B0",  "#DC050C", "#B2DF8A", "#B17BA6")

WholeTissueData_CD8T_ident$Status <- factor(WholeTissueData_CD8T_ident$Status,levels = c("Normal","AIP"))

tiff('SF7a.tiff',width = 1250,height = 700,res=150)
DimPlot(WholeTissueData_CD8T_ident,reduction = 'umap', label = F,repel = T,cols = .cluster_cols,group.by = "celltype2",split.by = "Status",pt.size = 0.5)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


PlotGenes <- c("AIF1","MAL","ACTN1","LDHB","NPM1","ARMH1","CCR7","LEF1","SELL","TCF7",
               "IFNG-AS1","INPP4B","FXYD2","CD84","CD44","CNN2","LYST","GPR183","GZMK","CD74",
               "ITGAE","CAPG","ITGA1","GPR25","SOCS1","JUN","FOS","CD9","ITM2C","CD69",
               "GNLY","FGFBP2","CX3CR1","FCGR3A","TBX21","KLRD1","ADGRG1","ITGB1","GZMB","GZMH",
               "CEBPD","TRAV1-2","NCR3","SLC4A10","ZBTB16","RORC","LTK","CXXC5","CCR6","KLRB1")

Sub_Data <- WholeTissueData_CD8T_ident

Sub_TF_Matrix <- Sub_Data@assays$RNA@data[PlotGenes,] %>% as.matrix()

Sub_TF_Matrix <- Sub_TF_Matrix[apply(Sub_TF_Matrix,1,sum) > 0,]

Sub_ClusterID <- lapply(split(Sub_Data@meta.data,list(Sub_Data$celltype2)), function(x)rownames(x))

Sub_TF_ClusterMean <- t(apply(Sub_TF_Matrix,1,function(x){
  lapply(Sub_ClusterID,function(ID){
    mean(x[ID])
  }) %>% unlist
}))

Sub_TF_ClusterMean <- Sub_TF_ClusterMean[apply(Sub_TF_ClusterMean,1,sum)>0,]
Sub_TF_ClusterMean[,1:ncol(Sub_TF_ClusterMean)] <- t(apply(Sub_TF_ClusterMean,1,scale))

Sub_barcolor <- .cluster_cols
names(Sub_barcolor) <- unique(Sub_Data$celltype2)
Sub_SampleBar <-  unique(Sub_Data$celltype2)
names(Sub_SampleBar) <-  unique(Sub_Data$celltype2)


library(ComplexHeatmap)
Sub_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Sub_SampleBar,
                                                col = list(bar = Sub_barcolor),show_legend = F,show_annotation_name = F)

pdf("SF7b.pdf",width = 8,height = 10)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),
                        show_row_names = F,show_column_names = T  ,column_names_rot = 45,top_annotation = Sub_ha_bar)
dev.off()



EachSampleClusterDisB <- lapply(unique(WholeTissueData_CD8T_ident@meta.data$orig.ident),function(i){
  seur <- subset(WholeTissueData_CD8T_ident,cells=rownames(WholeTissueData_CD8T_ident@meta.data[WholeTissueData_CD8T_ident@meta.data$orig.ident %in% i,]))
  seur@meta.data$celltype2 %>% table %>% data.frame %>% set_colnames(c("CellTypes","Number")) %>% dplyr::mutate(Per=100*Number/sum(Number))
})
names(EachSampleClusterDisB) <- unique(WholeTissueData_CD8T_ident@meta.data$orig.ident)
EachSampleClusterDisC <- dplyr::bind_rows(EachSampleClusterDisB) %>% dplyr::mutate(Samples=rep(names(EachSampleClusterDisB),times=unlist(lapply(EachSampleClusterDisB,nrow))))
EachSampleClusterDisC$Tissue <- gsub("N.","NI",EachSampleClusterDisC$Samples)
EachSampleClusterDisC$Tissue <- gsub("P.*","AIP",EachSampleClusterDisC$Tissue)
pdf('SF7c.pdf',width = 6.8,height = 4)
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

readr::write_rds(WholeTissueData_CD8T_ident,path="result/WholeTissueData_CD8T_ident",compress="gz")