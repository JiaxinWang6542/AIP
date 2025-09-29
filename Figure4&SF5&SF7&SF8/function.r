setwd("/work/wjx/scRNA-seq/AIP/figure/figure4")
library(Seurat)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)


WholeTissueData_CD4T_ident=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_CD4T_ident.rds.gz")
WholeTissueData_CD4T <- WholeTissueData_CD4T_ident


deg <- FindMarkers(WholeTissueData_CD4T,ident.1 = 'Tfh', only.pos = T)
deg_sig <- subset(deg, p_val_adj<0.05&avg_logFC>0.3)
ego_BP <- enrichGO(gene          = row.names(deg_sig),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 

ego_bp <- data.frame(ego_BP)
row.names(ego_bp)<- c(1:length(ego_bp$ID))
GODF <- ego_bp[c("10","22","28","30","32","41","66","68","85","103"),]
pdf("Figure4f.pdf",width = 12,height = 6)
ggplot(GODF,aes(x = -log10(p.adjust),y=Description,fill = Count))+
  geom_bar(stat = 'identity',width = 0.8)+
  scale_y_discrete(limits = rev(GODF$Description))+
  scale_fill_gradientn(colors =  colorRampPalette(c(RColorBrewer::brewer.pal(6,'Reds')),space='rgb')(10))+
  theme(axis.text=element_text(color = "black"),panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"),axis.title.y = element_blank(),axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12),axis.title.x = element_text(size=15))
dev.off()



####clusterProfiler GO 
WholeTissueData_CD8T_ident=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_CD8T_ident.rds.gz")
WholeTissueData_CD8T <- WholeTissueData_CD8T_ident
table(WholeTissueData_CD8T$celltype2)
table(WholeTissueData_CD8T$celltype)

deg <- FindMarkers(WholeTissueData_CD8T,ident.1 = 'Trm', only.pos = T)
deg_sig <- subset(deg, p_val_adj<0.05&avg_logFC>0.3)


ego_BP <- enrichGO(gene          = row.names(deg_sig),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 

ego_bp <- data.frame(ego_BP)
row.names(ego_bp)<- c(1:length(ego_bp$ID))
GODF <- ego_bp[c("5","10","12","13","15","20","29","58","64","65"),]
pdf("SF7d.pdf",width = 9,height = 6)
ggplot(GODF,aes(x = -log10(p.adjust),y=Description,fill = Count))+
  geom_bar(stat = 'identity',width = 0.8)+
  scale_y_discrete(limits = rev(GODF$Description))+
  scale_fill_gradientn(colors =  colorRampPalette(c(RColorBrewer::brewer.pal(6,'Reds')),space='rgb')(10))+
  theme(axis.text=element_text(color = "black"),panel.background = element_blank(),axis.line = element_line(size = 0.5, colour = "black"),axis.title.y = element_blank(),axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12),axis.title.x = element_text(size=15))
dev.off()