
setwd("/work/wjx/scRNA-seq/AIP/figure/figure5")
library(Seurat)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ggsci)
library(magrittr)
library(ComplexHeatmap)
#######TLS score
WholeTissueData <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure1/WholeTissueData_ident.rds.gz")
TLS_genes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11","CXCL13")

###addmodulescore
gene <- list(TLS_genes)
WholeTissueData <- AddModuleScore(WholeTissueData,gene,name="TLS_add")
WholeTissueData <- subset(WholeTissueData,TLS_add1<0.6)
ECMCompareDF <- WholeTissueData@meta.data[,c('TLS_add1','celltype','Status')]
ECMCompareDF$Status <- factor(ECMCompareDF$Status,levels=c("Normal","AIP"))
pdf('Figure5a.pdf',width =4,height =4)
ggplot(ECMCompareDF,aes(x =Status ,y = TLS_add1, fill = Status))+
  ylab("TLS Score")+
  geom_boxplot(stat = "boxplot",alpha=1,width=.6, outlier.shape = NA)+
  scale_fill_manual(values = c("#1f78b4","#CB181D"))+
  theme(panel.background = element_blank(),panel.grid = element_blank(),strip.background = element_blank(),strip.text = element_text(size = 12),
        axis.title.x = element_blank(),legend.position = 'none',axis.line.y.left = element_line(color = 'black'),axis.ticks.y.left = element_line(color = 'black'),
       axis.text.y.left = element_text(color = 'black'),axis.title.y.left = element_text(color = 'black'),axis.line.x = element_line(color = 'black'),axis.text.x = element_text(size=10,color = 'black'))+
  ggpubr::stat_compare_means(method = 't.test')
dev.off()





