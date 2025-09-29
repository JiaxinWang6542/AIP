setwd("/work/wjx/scRNA-seq/AIP/")
library(scRepertoire)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(patchwork)
library(ggsci)
library(ggalluvial)

DefineB <- readr::read_rds("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_T_ident.rds.gz")

PatientList <- list.files("/work/wjx/scRNA-seq/AIP/TCR/filter_contig_annotations") 

Mtcr_list <- lapply(PatientList,function(x){
  Xdata <- read.csv(file = paste("/work/wjx/scRNA-seq/AIP/TCR/filter_contig_annotations/",x,sep = ""), stringsAsFactors = F)
  return(Xdata)
})
names(Mtcr_list)=gsub("_filtered_contig_annotations.csv","",PatientList )
head(Mtcr_list[[1]])

barcode <- rownames(DefineB@meta.data)

combined_Mtcr <- combineTCR(Mtcr_list,samples = names(Mtcr_list),ID= names(Mtcr_list),removeNA=T,filterMulti=T)
names(combined_Mtcr) <- names(Mtcr_list)
head(combined_Mtcr[[1]])
combined_Mtcr <- lapply(names(combined_Mtcr),function(Name){
  x <- combined_Mtcr[[Name]]
  rx <- Mtcr_list[[Name]]
  x$barcode <- gsub("N[1-8]_N[1-8]_|P[1-9]_P[1-9]_|P1_1_P1_1_|P1_2_P1_2_|P10_P10_","",x$barcode)
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
combined_Mtcr <- addVariable(combined_Mtcr, variable.name = "patient", variables = c("N1","N2","N3","N4","N5","N6","N7","N8","P1","P1","P10","P2","P3","P4","P5","P6","P7","P8","P9"))
combined_Mtcr <- addVariable(combined_Mtcr, variable.name = "group", variables = c(rep("Normal",8),rep("AIP",11)))


data <- clonalQuant(combined_Mtcr, 
            cloneCall="aa", 
            chain = "both", 
            scale = TRUE,
            group.by = "patient",
            exportTable = T,
            palette = "viridis")
data$group <- gsub("N.*","Normal",gsub("P.*","AIP",data$patient))
pdf('SF8b.pdf',width = 4,height = 3.5)
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
                cloneCall = "aa",
                group.by="patient",
                metrics = c("shannon", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE"),
                exportTable = T)
data$group <- gsub("N.*","Normal",gsub("P.*","AIP",data$patient))
data <- data[,-8]
data_new <- gather(data,key="diversity",value="y",-"patient",-"group")

pdf('SF8a.pdf',width = 10,height = 7)
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
                                   cloneCall="aa", 
                             group.by = "sample", 
                                   proportion = F, 
                                   cloneSize=c(Small=1, Medium=5, Large=20))
table(DefineT$celltype1)
head(DefineT@meta.data)
tiff("SF8c.tiff",width = 1000,height = 700,res=150)
Seurat::DimPlot(DefineT, group.by = "cloneSize") +
  scale_color_manual(values=c("#DC050C","#FF7F00","#7BAFDE"))
dev.off()



pdf("SF8d.pdf",width = 8,height = 6)
clonalOccupy(DefineT, 
               x.axis = "celltype1", 
             proportion = T, 
             label = F)+
  scale_fill_manual(values=c("#DC050C","#FF7F00","#7BAFDE"))
dev.off()




pdf("SF8e.pdf",width = 8,height = 6)
clonalOverlap(DefineT, 
              cloneCall = "aa", 
              method = "cosine",
              group.by="celltype1",
              palette="Viridis")+
  scale_fill_gradientn(colors =  colorRampPalette(c(RColorBrewer::brewer.pal(8,'Reds')),space='rgb')(10),na.value = "white")
dev.off()

