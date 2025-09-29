
setwd("/work/wjx/scRNA-seq/AIP/figure/figure6")
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

#####CP DATA
#her CP 1-5
Xdata <- Read10X_h5("/work/wjx/scRNA-seq/AIP/CP_matrix/GSE165045_BL5_filtered_feature_bc_matrix.h5" )
cp_her5 <- CreateSeuratObject(Xdata[[1]],project="her5",min.cells=3,min.features = 200)
cp_her5

cp_her5$sample <- gsub(".+-","",rownames(cp_her5@meta.data))
#Control 6-8 
Xdata <- Read10X_h5("/work/wjx/scRNA-seq/AIP/CP_matrix/GSE165045_BL8_filtered_feature_bc_matrix.h5" )
cp_ctl3 <- CreateSeuratObject(Xdata[[1]],project="ctl3",min.cells=3,min.features = 200)
cp_ctl3

#idio CP 9-12
Xdata <- Read10X_h5("/work/wjx/scRNA-seq/AIP/CP_matrix/GSE165045_BL12_filtered_feature_bc_matrix.h5" )
cp_idio4 <- CreateSeuratObject(Xdata,project="idio4",min.cells=3,min.features = 200)
cp_idio4



WholeCPlist <- list(cp_her1,cp_her2,cp_her3,cp_her4,cp_her5,cp_ctl1,cp_ctl2,cp_ctl3,cp_idio1,cp_idio2,cp_idio3,cp_idio4)
names(WholeCPlist) <- c("her1","her2","her3","her4","her5","ctl1","ctl2","ctl3","idio1","idio2","idio3","idio4")
WholeTissueData_CP<- WholeCPlist[[1]]
WholeTissueData_CP<- merge(x = WholeTissueData_CP, y = WholeCPlist[2:length(WholeCPlist)],add.cell.ids=c(names(WholeCPlist)),merge.data=T,project = "SeuratProject")
head(rownames(WholeTissueData_CP@meta.data))
head(WholeTissueData_CP@meta.data)


WholeTissueData_CP@meta.data$scrublet_doublets <-   "False"
WholeTissueData_CP@meta.data$doubletdetection_doublets <-   "False"
WholeTissueData_CP@meta.data$Sample <-   WholeTissueData_CP@meta.data$orig.ident
table(WholeTissueData_CP@meta.data$Sample)
mito.genes <- grep(pattern = "^MT-", x = rownames(WholeTissueData_CP@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(WholeTissueData_CP@assays[["RNA"]][mito.genes, ])/Matrix::colSums(WholeTissueData_CP@assays[["RNA"]])
WholeTissueData_CP <- AddMetaData(object = WholeTissueData_CP, metadata = percent.mito, col.name = "percent.mito") 
Ribo.genes <- grep(pattern = "^RPL|^RPS", x = rownames(WholeTissueData_CP@assays[["RNA"]]), value = TRUE)
WholeTissueData_CP$percent.ribo <- Matrix::colSums(WholeTissueData_CP@assays[["RNA"]][Ribo.genes, ])/Matrix::colSums(WholeTissueData_CP@assays[["RNA"]])

HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
HB_m <- match(HB.genes_total,rownames(WholeTissueData_CP@assays$RNA))
HB.genes <- rownames(WholeTissueData_CP@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
WholeTissueData_CP[["percent.HB"]]<-PercentageFeatureSet(WholeTissueData_CP,features=HB.genes)

WholeTissueData_CP@meta.data$celltype <- "CPcell"
WholeTissueData_CP@meta.data$Status <- gsub("her.*","Her",WholeTissueData_CP@meta.data$orig.ident)
WholeTissueData_CP@meta.data$Status <- gsub("idio.*","Idio",WholeTissueData_CP@meta.data$Status)
WholeTissueData_CP@meta.data$Status <- gsub("ctl.*","Ctl",WholeTissueData_CP@meta.data$Status)
table(WholeTissueData_CP$Status)
WholeTissueData_CP@meta.data$celltype1 <- "CPcell"
head(WholeTissueData_CP@meta.data)
head(WholeTissueData_AIP@meta.data)
readr::write_rds(WholeTissueData_CP,path="WholeTissueData_CP_Raw.rds.gz",compress="gz")
rm(list=ls())
WholeTissueData_CP=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure6/WholeTissueData_CP_Raw.rds.gz")
head(WholeTissueData_CP@meta.data)

WholeTissueData_CP <- subset(WholeTissueData_CP, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mito >  -Inf & percent.mito < 0.2& percent.HB < 0.15 )



##pca-umap clustering
#wihtout harmony
table(WholeTissueData_CP$Status)
WholeTissueData_CP  <- NormalizeData(WholeTissueData_CP) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_CP)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_CP <- ScaleData(WholeTissueData_CP,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)

ElbowPlot(WholeTissueData_CP,ndims = 50)

WholeTissueData_CP <- FindNeighbors(WholeTissueData_CP,dims = 1:20)
WholeTissueData_CP  <- RunUMAP(object = WholeTissueData_CP, dims = 1:20) 
WholeTissueData_CP   <- FindClusters(WholeTissueData_CP , resolution = c(seq(0.1,0.8,0.1)), algorithm = 1)

DimPlot(WholeTissueData_CP,group.by = "orig.ident",label = T)
DimPlot(WholeTissueData_CP,group.by = "Status",label = T)


levels(Idents(WholeTissueData_CP))
Idents(WholeTissueData_CP) <- 'RNA_snn_res.0.4'


WholeTissueData_ident <- RenameIdents(WholeTissueData_CP, "0" = "T&NK","1" = "T&NK","2" = "T&NK","4" = "T&NK","12" = "T&NK","9"= "B","10"="Plasma","3"="Myeloid","6"="Myeloid","7"="Myeloid","8"="Myeloid","5"="Mast","15"="pDC","11"="Prolifer","16"="Ductal&Acinar","13"="Endothelial","14"="Fibroblast")

WholeTissueData_ident$celltype=WholeTissueData_ident@active.ident
WholeTissueData_ident_markers <- FindAllMarkers(WholeTissueData_ident,only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(WholeTissueData_ident_markers ,file = "bigcluster_markers_final.csv")

readr::write_rds(WholeTissueData_ident,path="./result/WholeTissueData_CP_ident.rds.gz",compress="gz")

WholeTissueData_ident <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure6/result/WholeTissueData_CP_ident.rds.gz")


WholeTissueData_main_ident <- WholeTissueData_ident
table(WholeTissueData_main_ident$celltype)

WholeTissueData_ident_AIP_main <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure1/WholeTissueData_ident.rds.gz")

colnames(WholeTissueData_ident_AIP_main@meta.data)
table(WholeTissueData_ident_AIP_main$celltype)
WholeTissueData_ident_AIP_main@meta.data <- WholeTissueData_ident_AIP_main@meta.data[,c(1,2,3,7,8,9,10,16,17)]
head(WholeTissueData_ident_AIP_main@meta.data)

colnames(WholeTissueData_main_ident@meta.data)
WholeTissueData_main_ident@meta.data <- WholeTissueData_main_ident@meta.data[,c(1,2,3,6,7,8,9,10,11)]
head(WholeTissueData_main_ident@meta.data)
WholeTissuedata_intergrated_main <- merge(x = WholeTissueData_ident_AIP_main, y = WholeTissueData_main_ident,merge.data=T,project = "SeuratProject")
readr::write_rds(WholeTissuedata_intergrated_main,path="/work/wjx/scRNA-seq/AIP/figure/figure6/result/WholeTissuedata_intergrated_main.rds.gz",compress="gz")


WholeTissuedata_intergrated_main  <- NormalizeData(WholeTissuedata_intergrated_main)
PlotGenes <- c("CCL20","CXCL2","CXCL3","CXCL8","CCL2","CCL3","CCL5",
               "CXCL9","CXCL10","CXCL11","CXCL13","CCL19","CCL21","IL21")
WholeTissuedata_intergrated_main$group <- gsub("Ctl|Normal","Normal",WholeTissuedata_intergrated_main$Status)
WholeTissuedata_intergrated_main$group <- gsub("Her|Idio","CP",WholeTissuedata_intergrated_main$group)
table(WholeTissuedata_intergrated_main$group,WholeTissuedata_intergrated_main$celltype)




WholeTissuedata_intergrated_main  <- NormalizeData(WholeTissuedata_intergrated_main) %>% FindVariableFeatures(nfeatures = 4000)
WholeTissuedata_intergrated_main <- ScaleData(WholeTissuedata_intergrated_main) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissuedata_intergrated_main,ndims = 50)
###harmony
table(WholeTissuedata_intergrated_main$orig.ident)
WholeTissuedata_intergrated_main <- WholeTissuedata_intergrated_main  %>% 
  RunHarmony("orig.ident",plot_convergence = TRUE)
WholeTissuedata_intergrated_main <- WholeTissuedata_intergrated_main %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)


table(WholeTissuedata_intergrated_main$celltype1)
WholeTissuedata_intergrated_main$celltype1 <- gsub("T&NK","NK",WholeTissuedata_intergrated_main$celltype)
WholeTissuedata_intergrated_main$celltype1 <- gsub("T|NK","T&NK",WholeTissuedata_intergrated_main$celltype1)
WholeTissuedata_intergrated_main$celltype1 <- gsub("DC|pDC","Myeloid",WholeTissuedata_intergrated_main$celltype1)
WholeTissuedata_intergrated_main$celltype1 <- gsub("Macrophage","Myeloid",WholeTissuedata_intergrated_main$celltype1)
WholeTissuedata_intergrated_main$celltype1 <- gsub("Neutrophil","Myeloid",WholeTissuedata_intergrated_main$celltype1)
WholeTissuedata_intergrated_main$celltype1 <- gsub("Monocyte","Myeloid",WholeTissuedata_intergrated_main$celltype1)
WholeTissuedata_intergrated_main$celltype1 <- gsub("Mast","Myeloid",WholeTissuedata_intergrated_main$celltype1)


WholeTissuedata_intergrated_main_subset <- subset(WholeTissuedata_intergrated_main,celltype1 %in% c("B","Myeloid","Plasma","Prolifer","T&NK"))
WholeTissuedata_intergrated_main_subset  <- NormalizeData(WholeTissuedata_intergrated_main_subset) %>% FindVariableFeatures(nfeatures = 4000)
WholeTissuedata_intergrated_main_subset <- ScaleData(WholeTissuedata_intergrated_main_subset) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissuedata_intergrated_main_subset,ndims = 50)


###harmony
table(WholeTissuedata_intergrated_main_subset$orig.ident)
WholeTissuedata_intergrated_main_subset <- WholeTissuedata_intergrated_main_subset  %>% 
  RunHarmony("orig.ident",plot_convergence = TRUE)
WholeTissuedata_intergrated_main_subset <- WholeTissuedata_intergrated_main_subset %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)


.cluster_cols <- c( "#DC050C","#1965B0","#33A02C","#FF7F00", "#882E72",
                    "#FB8072","#7BAFDE","#B2DF8A","#FDB462",  "#B17BA6", 
  "#E7298A", "#E78AC3","#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
  "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
  "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

pdf("SF9a.pdf",width = 6,height = 5)
DimPlot(WholeTissuedata_intergrated_main_subset,raster=FALSE,group.by = "celltype1",cols = .cluster_cols ,pt.size = 0.5)
dev.off()


pdf("SF9b.pdf",width = 6,height = 5)
DimPlot(WholeTissuedata_intergrated_main_subset,raster=FALSE,group.by = "group",pt.size = 0.5)
dev.off()


WholeTissuedata_intergrated_main_subset$celltype1 <- factor(WholeTissuedata_intergrated_main_subset$celltype1,levels=c("T&NK","B","Plasma","Myeloid","Prolifer"))
ClusterFreq <- WholeTissuedata_intergrated_main_subset@meta.data[,c("celltype1","group")] %>% table %>% 
  data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
ClusterPer$seurat_clusters=factor(ClusterPer$seurat_clusters,levels = c("T&NK","B","Plasma","Myeloid","Prolifer"))

pdf("SF9c.pdf",width = 7,height = 5)
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





Sub_Data <- WholeTissuedata_intergrated_main
Sub_Data$group <- factor(Sub_Data$group,levels = c("Normal","CP","AIP"))
DefaultAssay(Sub_Data) <- "RNA"

Sub_TF_Matrix <- Sub_Data@assays$RNA@data[PlotGenes,] %>% as.matrix()

Sub_TF_Matrix <- Sub_TF_Matrix[apply(Sub_TF_Matrix,1,sum) > 0,]

Sub_ClusterID <- lapply(split(Sub_Data@meta.data,list(Sub_Data$group)), function(x)rownames(x))

Sub_TF_ClusterMean <- t(apply(Sub_TF_Matrix,1,function(x){
  lapply(Sub_ClusterID,function(ID){
    mean(x[ID])
  }) %>% unlist
}))

Sub_TF_ClusterMean <- Sub_TF_ClusterMean[apply(Sub_TF_ClusterMean,1,sum)>0,]
Sub_TF_ClusterMean[,1:ncol(Sub_TF_ClusterMean)] <- t(apply(Sub_TF_ClusterMean,1,scale))

Sub_barcolor <- c( "#E7298A","#FF7F00","#DC050C")
names(Sub_barcolor) <- levels(Sub_Data$group)
Sub_SampleBar <-  levels(Sub_Data$group)
names(Sub_SampleBar) <-  levels(Sub_Data$group)


library(ComplexHeatmap)
Sub_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Sub_SampleBar,
                                                col = list(bar = Sub_barcolor),show_legend = F,show_annotation_name = F)

pdf("Figure6e.pdf",width = 4,height = 7)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu"))[2:9],space="rgb")(100),
                        show_row_names = T,show_column_names = T  ,column_names_rot = 45,top_annotation = Sub_ha_bar)
dev.off()





######T cell
table(WholeTissueData_ident$celltype)
WholeTissueData_T <- subset(WholeTissueData_ident,celltype %in% c("T&NK"))
table(WholeTissueData_T@meta.data$celltype1)

WholeTissueData_T  <- NormalizeData(WholeTissueData_T) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_T)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_T <- ScaleData(WholeTissueData_T,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissueData_T,ndims = 50)

WholeTissueData_harmony_T <- WholeTissueData_T %>% 
RunHarmony("Status", plot_convergence = TRUE)

WholeTissueData_harmony_T <- WholeTissueData_harmony_T %>% 
RunUMAP(reduction = "harmony", dims = 1:20) %>% 
FindNeighbors(reduction = "harmony", dims = 1:20)
WholeTissueData_harmony_T  <- FindClusters(WholeTissueData_harmony_T , resolution = c(seq(0.1,1,0.1)), algorithm = 1)

clustree(WholeTissueData_harmony_T @meta.data,prefix = "RNA_snn_res.")


Idents(WholeTissueData_harmony_T) <- "RNA_snn_res.0.8"



WholeTissueData_T_ident <- RenameIdents(WholeTissueData_harmony_T, "2" = "CD8+ T","3" = "CD8+ T","4" = "CD8+ T","8" = "CD8+ T","13" = "CD8+ T","14" = "CD8+ T","15" = "CD8+ T",
                                        "0" = "CD4+ T","1" = "CD4+ T","5" = "CD4+ T","6" = "CD4+ T","7" = "CD4+ T","9" = "CD4+ T","12" = "CD4+ T",
                                        "10" = "other","11" = "other","16" = "other")

DimPlot(WholeTissueData_T_ident)
WholeTissueData_T_ident$celltype <- WholeTissueData_T_ident@active.ident
readr::write_rds(WholeTissueData_T_ident,path="./result/WholeTissueData_CP_T_ident.rds.gz",compress="gz")


######CD4T
WholeTissueData_CD4T_ident <- subset(WholeTissueData_T_ident,celltype=="CD4+ T")
table(WholeTissueData_CD4T_ident$celltype)
WholeTissueData_ident_AIP_CD4T <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_CD4T_ident.rds.gz")

WholeTissueData_ident_AIP_CD4T$celltype1 <- WholeTissueData_ident_AIP_CD4T$celltype2
WholeTissueData_ident_AIP_CD4T@meta.data <- WholeTissueData_ident_AIP_CD4T@meta.data[,c(1,2,3,7,8,9,10,16,17,29)]
head(WholeTissueData_ident_AIP_CD4T@meta.data)

colnames(WholeTissueData_CD4T_ident@meta.data)
WholeTissueData_CD4T_ident@meta.data <- WholeTissueData_CD4T_ident@meta.data[,c(1,2,3,6,7,8,9,10,11,12)]
head(WholeTissueData_CD4T_ident@meta.data)

c1 <- which(rownames(WholeTissueData_ident_AIP_CD4T) %in% rownames(WholeTissueData_CD4T_ident))
WholeTissueData_ident_AIP_CD4T <- subset(WholeTissueData_ident_AIP_CD4T,features=rownames(WholeTissueData_ident_AIP_CD4T)[c1])
c1 <- which(rownames(WholeTissueData_CD4T_ident) %in% rownames(WholeTissueData_ident_AIP_CD4T))
WholeTissueData_CD4T_ident<- subset(WholeTissueData_CD4T_ident,features=rownames(WholeTissueData_CD4T_ident)[c1])


WholeTissueData_list_CD4T <- list(WholeTissueData_ident_AIP_CD4T,WholeTissueData_CD4T_ident)
names(WholeTissueData_list_CD4T) <- c('AIP','CP')



#####transfer anchors#####

WholeTissueData_list_CD4T <- lapply(WholeTissueData_list_CD4T, function(x){
  x<- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 4000)
  return(x)
})


#Identify anchors(AIP reference)
transfer.anchors <- FindTransferAnchors(reference = WholeTissueData_list_CD4T[[1]],query = WholeTissueData_list_CD4T[[2]],
                                        reference.assay = 'RNA',query.assay = 'RNA',reduction = 'pcaproject',reference.reduction="pca")

#label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors,refdata = WholeTissueData_list_CD4T[[1]]$celltype1,
                                     dims = 1:30)

WholeTissueData_list_CD4T[[2]] <- AddMetaData(WholeTissueData_list_CD4T[[2]],metadata = celltype.predictions)
table(WholeTissueData_list_CD4T[[2]]$Status)
WholeTissueData_list_CD4T[[2]]$celltype1 <- WholeTissueData_list_CD4T[[2]]$predicted.id




####merge
WholeTissuedata_intergrated_CD4T <- merge(x = WholeTissueData_list_CD4T[[1]], y = WholeTissueData_list_CD4T[[2]],merge.data=T,project = "SeuratProject")
table(WholeTissuedata_intergrated_CD4T$celltype1,WholeTissuedata_intergrated_CD4T$Status)
WholeTissuedata_intergrated_CD4T$group <- gsub("Her|Idio","CP",WholeTissuedata_intergrated_CD4T$Status)
WholeTissuedata_intergrated_CD4T$group <- gsub("Ctl|Normal","Normal",WholeTissuedata_intergrated_CD4T$group)

readr::write_rds(WholeTissuedata_intergrated_CD4T,path="./result/WholeTissuedata_intergrated_CD4T.rds.gz",compress="gz")


#######ROE
cluster.table <- table(WholeTissuedata_intergrated_CD4T$celltype1, WholeTissuedata_intergrated_CD4T$group)
cluster.table %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = purrr::pmap_dbl(
      list(
        .x = Var1,
        .y = Var2,
        .z = Freq
      ),
      .f = function(.x, .y, .z){
        a <- .z
        b <- sum(cluster.table[,.y]) - a
        c <- sum(cluster.table[.x,]) - a
        d <- sum(cluster.table) - a - b - c
        
        #o <- fisher.test(matrix(c(a, b,c, d), ncol = 2), alternative = "greater")
        #o$estimate
        o <- chisq.test(matrix(c(a, b, c, d), ncol = 2))
        oe <- o$observed/o$expected
        oe[1,1]
      }
    )
  ) -> enrich.res
enrich.res$p.value <- round(enrich.res$p.value,2)

pdf("Figure6b.pdf",width = 4.5,height = 5)
enrich.res %>%
  dplyr::rename(`Ro/e` = p.value) %>%
  ggplot(aes( Var2, Var1,fill = `Ro/e`)) +
  geom_tile(colour = "white", lwd = 0.8) +
  geom_text(aes(label=`Ro/e`),color='black') +
  scale_y_discrete(limits=rev(c("Tn","LIMS1+ Tm","TIMP1+ Tm","Tem","Temra","TNFRSF18- Treg","TNFRSF18+ Treg","Tfh","ISG+ T")))+
  scale_fill_continuous(limits=c(0,2.5))+
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black" ,angle = 0)) +
  scale_fill_distiller(palette = "Spectral") 
dev.off()


###CCA
WholeTissuedata_CD4T.list <- SplitObject(WholeTissuedata_intergrated_CD4T,split.by = "Status")
for (i in 1:length(WholeTissuedata_CD4T.list)) {
  WholeTissuedata_CD4T.list[[i]] <- NormalizeData(WholeTissuedata_CD4T.list[[i]],verbose = F)
  WholeTissuedata_CD4T.list[[i]] <- FindVariableFeatures(WholeTissuedata_CD4T.list[[i]],selection.method = 'vst',nfeatures = 4000,verbose = F)
}

CD4T.anchors <- FindIntegrationAnchors(object.list = WholeTissuedata_CD4T.list,dims = 1:50)
CD4T.integrated <- IntegrateData(anchorset = CD4T.anchors,dims = 1:50)

DefaultAssay(CD4T.integrated) <- 'integrated'
CD4T.integrated <- ScaleData(CD4T.integrated,verbose = F)
CD4T.integrated <- RunPCA(CD4T.integrated,npcs = 50,verbose = F)
CD4T.integrated <- RunUMAP(CD4T.integrated,reduction = 'pca',dims = 1:20)
.cluster_cols <- c(
  "#FB8072", "#1965B0", "#7BAFDE",  "#33A02C", "#B2DF8A","#FF7F00", "#FDB462","#DC050C","#882E72")
pdf("SF9f.pdf",width = 4.5,height = 5)
DimPlot(CD4T.integrated,split.by = "Status",label = F,cols = .cluster_cols)
dev.off()




######B cell
table(WholeTissueData_ident$celltype)
WholeTissueData_B_ident <- subset(WholeTissueData_ident,celltype %in% c("B","Plasma"))
table(WholeTissueData_B_ident$celltype)
WholeTissueData_ident_AIP_B <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")

colnames(WholeTissueData_ident_AIP_B@meta.data)
table(WholeTissueData_ident_AIP_B$celltype1)
WholeTissueData_ident_AIP_B@meta.data <- WholeTissueData_ident_AIP_B@meta.data[,c(1,2,3,7,8,9,10,16,17,32)]
head(WholeTissueData_ident_AIP_B@meta.data)

colnames(WholeTissueData_B_ident@meta.data)
WholeTissueData_B_ident@meta.data <- WholeTissueData_B_ident@meta.data[,c(1,2,3,6,7,8,9,10,11,12)]
head(WholeTissueData_B_ident@meta.data)

c1 <- which(rownames(WholeTissueData_ident_AIP_B) %in% rownames(WholeTissueData_B_ident))
WholeTissueData_ident_AIP_B <- subset(WholeTissueData_ident_AIP_B,features=rownames(WholeTissueData_ident_AIP_B)[c1])
c1 <- which(rownames(WholeTissueData_B_ident) %in% rownames(WholeTissueData_ident_AIP_B))
WholeTissueData_B_ident<- subset(WholeTissueData_B_ident,features=rownames(WholeTissueData_B_ident)[c1])


WholeTissueData_list_B <- list(WholeTissueData_ident_AIP_B,WholeTissueData_B_ident)
names(WholeTissueData_list_B) <- c('AIP','CP')



#####transfer anchors#####
WholeTissueData_list_B <- lapply(WholeTissueData_list_B, function(x){
  x<- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 4000)
  return(x)
})


#Identify anchors(AIP reference)
transfer.anchors <- FindTransferAnchors(reference = WholeTissueData_list_B[[1]],query = WholeTissueData_list_B[[2]],
                                        reference.assay = 'RNA',query.assay = 'RNA',reduction = 'pcaproject',reference.reduction="pca")

#label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors,refdata = WholeTissueData_list_B[[1]]$celltype1,
                                     dims = 1:30)

WholeTissueData_list_B[[2]] <- AddMetaData(WholeTissueData_list_B[[2]],metadata = celltype.predictions)
table(WholeTissueData_list_B[[2]]$Status)
WholeTissueData_list_B[[2]]$celltype1 <- WholeTissueData_list_B[[2]]$predicted.id




####merge
WholeTissuedata_intergrated_B <- merge(x = WholeTissueData_list_B[[1]], y = WholeTissueData_list_B[[2]],merge.data=T,project = "SeuratProject")
table(WholeTissuedata_intergrated_B$celltype1,WholeTissuedata_intergrated_B$Status)



#######ROE
cluster.table <- table(WholeTissuedata_intergrated_B$celltype1, WholeTissuedata_intergrated_B$group)
cluster.table %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = purrr::pmap_dbl(
      list(
        .x = Var1,
        .y = Var2,
        .z = Freq
      ),
      .f = function(.x, .y, .z){
        a <- .z
        b <- sum(cluster.table[,.y]) - a
        c <- sum(cluster.table[.x,]) - a
        d <- sum(cluster.table) - a - b - c
        
        #o <- fisher.test(matrix(c(a, b,c, d), ncol = 2), alternative = "greater")
        #o$estimate
        o <- chisq.test(matrix(c(a, b, c, d), ncol = 2))
        oe <- o$observed/o$expected
        oe[1,1]
      }
    )
  ) -> enrich.res
enrich.res$p.value <- round(enrich.res$p.value,2)

pdf("Figure6a.pdf",width = 4,height = 4)
enrich.res %>%
  dplyr::rename(`Ro/e` = p.value) %>%
  ggplot(aes( Var2, Var1,fill = `Ro/e`)) +
  geom_tile(colour = "white", lwd = 0.8) +
  geom_text(aes(label=`Ro/e`),color='black') +
  scale_y_discrete(limits=rev(c("Naive B","Memory B","ABC(IgD+)","ABC(IgD-)","GC B","Plasma A","Plasma G")))+
  scale_fill_continuous(limits=c(0,2.5))+
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black" ,angle = 0)) +
  scale_fill_distiller(palette = "Spectral") 
dev.off()


###CCA
WholeTissuedata_B.list <- SplitObject(WholeTissuedata_intergrated_B,split.by = "Status")
for (i in 1:length(WholeTissuedata_B.list)) {
  WholeTissuedata_B.list[[i]] <- NormalizeData(WholeTissuedata_B.list[[i]],verbose = F)
  WholeTissuedata_B.list[[i]] <- FindVariableFeatures(WholeTissuedata_B.list[[i]],selection.method = 'vst',nfeatures = 4000,verbose = F)
}

B.anchors <- FindIntegrationAnchors(object.list = WholeTissuedata_B.list,dims = 1:50)
B.integrated <- IntegrateData(anchorset = B.anchors,dims = 1:50)

DefaultAssay(B.integrated) <- 'integrated'
B.integrated <- ScaleData(B.integrated,verbose = F)
B.integrated <- RunPCA(B.integrated,npcs = 50,verbose = F)
B.integrated <- RunUMAP(B.integrated,reduction = 'pca',dims = 1:20)
.cluster_cols <- c("#B2DF8A","#1965B0","#FB8072", "#DC050C","#FF7F00","#882E72","#B17BA6")
pdf("SF9d.pdf",width = 4.5,height = 5)
DimPlot(B.integrated,split.by = "Status",label = F,cols = .cluster_cols)
dev.off()





#########myeloid
table(WholeTissueData_ident$celltype)
WholeTissueData_myeloid <- subset(WholeTissueData_ident,celltype %in% c("Myeloid"))
table(WholeTissueData_myeloid@meta.data$celltype1)
WholeTissueData_myeloid  <- NormalizeData(WholeTissueData_myeloid) %>% FindVariableFeatures(nfeatures = 4000)
VGENES=VariableFeatures(WholeTissueData_myeloid)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HBB|^HBA",VGENES)])
WholeTissueData_myeloid <- ScaleData(WholeTissueData_myeloid,features =VGENES,vars.to.regress =c("percent.mito","nFeature_RNA","nCount_RNA","percent.ribo","percent.HB")) %>% RunPCA(npcs = 50, verbose = FALSE)
ElbowPlot(WholeTissueData_myeloid,ndims = 50)

WholeTissueData_harmony_myeloid <- WholeTissueData_myeloid %>% 
  RunHarmony("Status", plot_convergence = TRUE)

WholeTissueData_harmony_myeloid <- WholeTissueData_harmony_myeloid %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
WholeTissueData_harmony_myeloid  <- FindClusters(WholeTissueData_harmony_myeloid , resolution = c(seq(0.1,1,0.1)), algorithm = 1)

clustree(WholeTissueData_harmony_myeloid@meta.data,prefix = "RNA_snn_res.")


Idents(WholeTissueData_harmony_myeloid) <- "RNA_snn_res.0.7"
markers <- FindAllMarkers(WholeTissueData_harmony_myeloid,only.pos = T,min.pct = 0.25)


WholeTissueData_myeloid_ident <- RenameIdents(WholeTissueData_harmony_myeloid, "1" = "macrophage","4" = "macrophage","5" = "macrophage","6" = "macrophage","7" = "macrophage","8" = "macrophage","9" = "macrophage",
                                              "2" = "monocyte","10" = "monocyte","0" = "DC","11" = "other","3" = "other"
                                        )

WholeTissueData_myeloid_ident$celltype <- WholeTissueData_myeloid_ident@active.ident


######macrophage
WholeTissueData_macrophage_ident <- subset(WholeTissueData_myeloid_ident,celltype=="macrophage")
table(WholeTissueData_macrophage_ident$celltype)
WholeTissueData_ident_AIP_macrophage <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure3/WholeTissueData_macrophage_ident.rds.gz")

colnames(WholeTissueData_ident_AIP_macrophage@meta.data)
table(WholeTissueData_ident_AIP_macrophage$celltype1)
WholeTissueData_ident_AIP_macrophage@meta.data <- WholeTissueData_ident_AIP_macrophage@meta.data[,c(1,2,3,7,8,9,10,16,17,29)]
head(WholeTissueData_ident_AIP_macrophage@meta.data)

colnames(WholeTissueData_macrophage_ident@meta.data)
WholeTissueData_macrophage_ident@meta.data <- WholeTissueData_macrophage_ident@meta.data[,c(1,2,3,6,7,8,9,10,11,12)]
head(WholeTissueData_macrophage_ident@meta.data)

c1 <- which(rownames(WholeTissueData_ident_AIP_macrophage) %in% rownames(WholeTissueData_macrophage_ident))
WholeTissueData_ident_AIP_macrophage <- subset(WholeTissueData_ident_AIP_macrophage,features=rownames(WholeTissueData_ident_AIP_macrophage)[c1])
c1 <- which(rownames(WholeTissueData_macrophage_ident) %in% rownames(WholeTissueData_ident_AIP_macrophage))
WholeTissueData_macrophage_ident<- subset(WholeTissueData_macrophage_ident,features=rownames(WholeTissueData_macrophage_ident)[c1])


WholeTissueData_list_macrophage <- list(WholeTissueData_ident_AIP_macrophage,WholeTissueData_macrophage_ident)
names(WholeTissueData_list_macrophage) <- c('AIP','CP')



#####transfer anchors#####
#Preprocess accroding to paper
WholeTissueData_list_macrophage <- lapply(WholeTissueData_list_macrophage, function(x){
  x<- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 4000)
  return(x)
})


#Identify anchors(AIP reference)
transfer.anchors <- FindTransferAnchors(reference = WholeTissueData_list_macrophage[[1]],query = WholeTissueData_list_macrophage[[2]],
                                        reference.assay = 'RNA',query.assay = 'RNA',reduction = 'pcaproject',reference.reduction="pca")

#label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors,refdata = WholeTissueData_list_macrophage[[1]]$celltype1,
                                     dims = 1:30)

WholeTissueData_list_macrophage[[2]] <- AddMetaData(WholeTissueData_list_macrophage[[2]],metadata = celltype.predictions)
table(WholeTissueData_list_macrophage[[2]]$Status)
WholeTissueData_list_macrophage[[2]]$celltype1 <- WholeTissueData_list_macrophage[[2]]$predicted.id




####merge
WholeTissuedata_intergrated_macrophage <- merge(x = WholeTissueData_list_macrophage[[1]], y = WholeTissueData_list_macrophage[[2]],merge.data=T,project = "SeuratProject")
table(WholeTissuedata_intergrated_macrophage$celltype1,WholeTissuedata_intergrated_macrophage$Status)
WholeTissuedata_intergrated_macrophage$group <- gsub("Ctl|Normal","Normal",WholeTissuedata_intergrated_macrophage$Status)
WholeTissuedata_intergrated_macrophage$group <- gsub("Her|Idio","CP",WholeTissuedata_intergrated_macrophage$group)
table(WholeTissuedata_intergrated_macrophage$group)
readr::write_rds(WholeTissuedata_intergrated_macrophage,path="./result/WholeTissuedata_intergrated_macrophage_20240810.rds.gz",compress="gz")

#######ROE
cluster.table <- table(WholeTissuedata_intergrated_macrophage$celltype1, WholeTissuedata_intergrated_macrophage$group)
cluster.table %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = purrr::pmap_dbl(
      list(
        .x = Var1,
        .y = Var2,
        .z = Freq
      ),
      .f = function(.x, .y, .z){
        a <- .z
        b <- sum(cluster.table[,.y]) - a
        c <- sum(cluster.table[.x,]) - a
        d <- sum(cluster.table) - a - b - c
        
        #o <- fisher.test(matrix(c(a, b,c, d), ncol = 2), alternative = "greater")
        #o$estimate
        o <- chisq.test(matrix(c(a, b, c, d), ncol = 2))
        oe <- o$observed/o$expected
        oe[1,1]
      }
    )
  ) -> enrich.res
enrich.res$p.value <- round(enrich.res$p.value,2)
pdf("Figure6c.pdf",width = 4.8,height = 3)
enrich.res %>%
  dplyr::rename(`Ro/e` = p.value) %>%
  ggplot(aes( Var2, Var1,fill = `Ro/e`)) +
  geom_tile(colour = "white", lwd = 0.8) +
  geom_text(aes(label=`Ro/e`),color='black') +
  scale_y_discrete(limits=rev(c("CXCL13+ Macrophage","CXCL9+ Macrophage","TREM2+ Macrophage","FOLR2+ Macrophage","FCN1+ Macrophage")))+
  scale_fill_continuous(limits=c(0,2.5))+
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black" ,angle = 0)) +
  scale_fill_distiller(palette = "Spectral") 
dev.off()


###CCA
WholeTissuedata_Macro.list <- SplitObject(WholeTissuedata_intergrated_macrophage,split.by = "Status")
for (i in 1:length(WholeTissuedata_Macro.list)) {
  WholeTissuedata_Macro.list[[i]] <- NormalizeData(WholeTissuedata_Macro.list[[i]],verbose = F)
  WholeTissuedata_Macro.list[[i]] <- FindVariableFeatures(WholeTissuedata_Macro.list[[i]],selection.method = 'vst',nfeatures = 4000,verbose = F)
}

Macro.anchors <- FindIntegrationAnchors(object.list = WholeTissuedata_Macro.list,dims = 1:50)
Macro.integrated <- IntegrateData(anchorset = Macro.anchors,dims = 1:50)

DefaultAssay(Macro.integrated) <- 'integrated'
Macro.integrated <- ScaleData(Macro.integrated,verbose = F)
Macro.integrated <- RunPCA(Macro.integrated,npcs = 50,verbose = F)
Macro.integrated <- RunUMAP(Macro.integrated,reduction = 'pca',dims = 1:20)
.cluster_cols <- c("#DC050C","#FF7F00","#33A02C","#1965B0","#B17BA6")
pdf("SF9e.pdf",width = 4.5,height = 5)
DimPlot(Macro.integrated,split.by = "Status",label = F,cols = .cluster_cols)
dev.off()

