setwd("/work/wjx/scRNA-seq/AIP/figure/figure2")
library(Seurat)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)


WholeTissueData_B_ident=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")

deg <- FindMarkers(WholeTissueData_B,ident.1 = 'ABC(IgD-)', only.pos = T)
deg_sig <- subset(deg, p_val_adj<0.1&avg_logFC>0.3)


ego_BP <- enrichGO(gene          = row.names(deg_sig),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 

ego_bp <- data.frame(ego_BP)
row.names(ego_bp)<- c(1:length(ego_bp$ID))
GODF <- ego_bp[c("15","32","47","48","56","98","106","189","198","286"),]
pdf("Figure2e.pdf",width = 12,height = 6)
ggplot(GODF,aes(x = -log10(p.adjust),y=Description,fill = Count))+
  geom_bar(stat = 'identity',width = 0.8)+
  scale_y_discrete(limits = rev(GODF$Description))+
  scale_fill_gradientn(colors =  colorRampPalette(c(RColorBrewer::brewer.pal(6,'Reds')),space='rgb')(10))+
  theme(axis.text=element_text(color = "black"),panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"),axis.title.y = element_blank(),axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12),axis.title.x = element_text(size=15))
dev.off()



PlotGenes <- c("CR2","CXCR5","BCL6","ADA",	"ADAM17","FOXJ1",	"KLHL6","TBX21","ITGAX","IGHD","IGHM","IGHA1","IGHG1","IGHG2","IGHG3","IGHE","IGHA2","IGHG4","AICDA","CXCR3","CXCR4")
Sub_Data <- subset(WholeTissueData_B_ident,celltype1 %in% c("GC B","ABC(IgD+)","ABC(IgD-)"))
Sub_Data$celltype1 <- factor(Sub_Data$celltype1,levels = c("GC B","ABC(IgD+)","ABC(IgD-)"))
DefaultAssay(Sub_Data) <- "RNA"

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

Sub_barcolor <- c( "#E7298A","#FF7F00","#DC050C")
names(Sub_barcolor) <- unique(Sub_Data$celltype1)
Sub_SampleBar <-  unique(Sub_Data$celltype1)
names(Sub_SampleBar) <-  unique(Sub_Data$celltype1)


library(ComplexHeatmap)
Sub_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Sub_SampleBar,
                                                col = list(bar = Sub_barcolor),show_legend = F,show_annotation_name = F)

pdf("Figure2f.pdf",width = 3.3,height = 7)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),
                        show_row_names = T,show_column_names = T  ,column_names_rot = 45,top_annotation = Sub_ha_bar)
dev.off()


