setwd("/work/wjx/scRNA-seq/AIP/figure/figure4")
library(Seurat)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ggsci)
library(magrittr)
library(ComplexHeatmap)
WholeTissueData_B <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")
WholeTissueData_B$celltype2 <- WholeTissueData_B$celltype1

WholeTissueData_CD4T <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_CD4T_ident.rds.gz")

a <- list(WholeTissueData_B,WholeTissueData_CD4T)
names(a) <- c("Bcell", "CD4T")

freq <- lapply(a,function(x){
  EachSampleClusterDisB <-lapply(unique(x@meta.data$orig.ident),function(i){
    seur <- subset(x,cells=rownames(x@meta.data[x@meta.data$orig.ident %in% i,]))
    seur@meta.data$celltype2 %>% table %>% data.frame %>% set_colnames(c("CellTypes","Number")) %>% dplyr::mutate(Per=100*Number/sum(Number))
  })
  names(EachSampleClusterDisB) <- unique(x@meta.data$orig.ident)
  EachSampleClusterDisC <- dplyr::bind_rows(EachSampleClusterDisB) %>% dplyr::mutate(Samples=rep(names(EachSampleClusterDisB),times=unlist(lapply(EachSampleClusterDisB,nrow))))
  EachSampleClusterDisC$Tissue <- gsub("N.*","Normal",EachSampleClusterDisC$Samples)
  EachSampleClusterDisC$Tissue <- gsub("P.*","AIP",EachSampleClusterDisC$Tissue)
  return(EachSampleClusterDisC)})

freq_data <-bind_rows(freq)
freq_data <- freq_data[,c(1,3,4)]
data <- spread(freq_data,key="CellTypes",value="Per",fill=0)
rownames(data) <- data[,1]
data <- data[,-1]
data <- data[,c(4,8:16)]
expr <- cor(data,method = "spearman")

library(ComplexHeatmap)
pdf("Figure4e.pdf",width = 5,height = 4.7)
Heatmap(expr,
        cluster_columns = T,clustering_distance_columns="euclidean",
        cluster_rows = T,clustering_distance_rows="euclidean",
        show_row_dend = F,column_dend_height = unit(20, "mm"),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = NA),
        show_column_names = T,
        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")),space="rgb")(100),
        heatmap_legend_param = list(title = "Correlation",title_gp = gpar(fontsize = 12, fontface = "bold"),legend_height=unit(50, "mm")),
)
dev.off()





