setwd("/work/wjx/scRNA-seq/AIP/figure/figure3")
library(Seurat)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr) 
library(fgsea)
library(future)
WholeTissueData_myeloid_ident=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure3/WholeTissueData_macrophage_ident.rds.gz")


#GSEA
plan(multisession,workers=40)
options(future.globals.maxSize= 3145728000)

mdb <- msigdbr(species = "Homo sapiens", category = "C5")
table(mdb$gs_subcat)
mdb <- mdb[grep("^GO:BP",mdb$gs_subcat),]
fgsea_sets<- mdb %>% split(x = .$gene_symbol, f = .$gs_name)
table(WholeTissueData_myeloid_ident$celltype1)
Idents(WholeTissueData_myeloid_ident) <- "Status"
cmarkers <- FindMarkers(WholeTissueData_myeloid_ident,ident.1 = "AIP",ident.2 = "Normal",min.pct = 0.01,logfc.threshold = 0.01)

head(cmarkers)
cmarkers$genes = rownames(cmarkers)
cluster0.genes<- cmarkers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
ranks<- deframe(cluster0.genes)
head(ranks)

geneList= cmarkers$avg_log2FC
names(geneList)= rownames(cmarkers)
geneList=sort(geneList,decreasing = T)
head(geneList)

m_df <- msigdbr(species = "Homo sapiens",category ="C5" )
mf_df <- mdb %>% dplyr::select(gs_name,gene_symbol)
colnames(mf_df) <- c("term","gene")

egmt <- GSEA(ranks,TERM2GENE=mf_df,pvalueCutoff = 0.05,pAdjustMethod = "fdr")
sortegmt<- egmt[order(egmt$NES, decreasing = T),]
sortegmt<- sortegmt[sortegmt$p.adjust <0.05,]
row.names(sortegmt) <- c(1:length(sortegmt$ID))
geneset_plot <- c(paste("GOBP_",toupper("leukocyte_migration"),sep = ""),paste("GOBP_",toupper("leukocyte_chemotaxis"),sep = ""))
library(enrichplot)
pdf("result/Fig3d.pdf",width = 7,height =4)
gseaplot2(egmt,geneSetID = geneset_plot,pvalue_table=T)
dev.off()






