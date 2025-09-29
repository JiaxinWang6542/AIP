setwd("/work/wjx/scRNA-seq/AIP/figure/figure2")
library(Seurat)
library(ggplot2)
library(pagoda2)
library(SeuratWrappers)
library(monocle)
library(htmlwidgets)
library(DDRTree)
library(pheatmap)
library(harmony)
library(tidyverse)

options(stringsAsFactors = F)
.cluster_cols <- c("#B2DF8A","#1965B0","#FB8072", "#DC050C","#FF7F00","#882E72","#B17BA6")

set.seed(123)

WholeTissueData_B <- readr::read_rds("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")

####Monocle Analysis
WholeTissueData_ABC <- subset(WholeTissueData_B,celltype1 %in% c("ABC(IgD+)","ABC(IgD-)","Plasma G")&Status=="AIP")
table(WholeTissueData_ABC$Status)
data <- as(as.matrix(WholeTissueData_ABC@assays$RNA@counts),'dgCMatrix')
pd <- new('AnnotatedDataFrame',data=WholeTissueData_ABC@meta.data)
fData <- data.frame(gene_short_name=row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame',data=fData)


HSMM <- monocle::newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = VGAM::negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)


HSMM <- detectGenes(HSMM,min_expr = 0.1)
print(head(fData(HSMM)))
expresed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
print(head(pData(HSMM)))


pie <- ggplot(pData(HSMM),aes(x=factor(1),fill=factor(celltype1)))+geom_bar(width = 1)
pie+coord_polar(theta = 'y')+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())

diff_test_res <- differentialGeneTest(HSMM[expresed_genes,],fullModelFormulaStr = '~celltype1',cores=30)
deg <- subset(diff_test_res,qval < 0.01)
deg <- deg[order(deg$qval,decreasing = F),]
head(deg)
ordering_genes <- row.names(deg)

ordering_genes <- row.names(deg)[1:2000]
HSMM_B <- setOrderingFilter(HSMM,ordering_genes)
plot_ordering_genes(HSMM_B)


HSMM_B <- reduceDimension(HSMM_B,max_components=2,method='DDRTree')

monocle::plot_cell_trajectory(HSMM_B,color_by = 'celltype1')+
     scale_color_manual(values = c("#FF7F00", "#DC050C","#1965B0"),name="")+
     theme(legend.position = 'right')
 


HSMM_B <- monocle::orderCells(HSMM_B,reverse = T)

plot_cell_trajectory(cds = HSMM_B,color_by = 'Pseudotime')+
  scale_color_gradient(low = '#D3D3D3',high = 'PURPLE')+
  theme(legend.position = 'right')

plot_cell_trajectory(HSMM_B,color_by = 'State')+
  theme(legend.position = 'right')

HSMM_B <- monocle::orderCells(HSMM_B,root_state = "5")




#plot
pdf('Figure2g_1.pdf',width = 5.5,height = 4)
monocle::plot_cell_trajectory(cds = HSMM_B,color_by = 'Pseudotime')+
  scale_color_gradient(low = '#D3D3D3',high = 'PURPLE')+
  theme(legend.position = 'right')
dev.off()

pdf('Figure2g_2.pdf',width = 6,height = 5.5)
monocle::plot_cell_trajectory(HSMM_B,color_by = 'celltype1')+
  scale_color_manual(values = c("#FF7F00", "#DC050C","#1965B0"),name="")+
  facet_wrap("~celltype1",nrow = 2)+
  theme(legend.position = 'none')
dev.off()


diff_test_res <- differentialGeneTest(HSMM_B[ordering_genes,],fullModelFormulaStr = '~sm.ns(Pseudotime)',cores=30)
select_genes <- c("TBX21","LRRK1","BATF","ITGAX","LMO2","RHOC","ID3","DHRS9","IL13RA1","IER5",
                  "CD53","SELL","CCR1","S100A4","MS4A1","AICDA","CXCR4","PFN1","PRDX1","SERF2",
                  "IGHE","IGHM","CD38","SSR3","LMAN1","HM13","RPN2","MANF","PDXK","CRELD2",
                 "MZB1","XBP1","JCHAIN","IGHG1","IGHG4","IGLC2","IGKC","IGHG2","IGHG3","SDC1")


p <- monocle::plot_pseudotime_heatmap(HSMM_B[select_genes,],
                                      hmcols = colorRampPalette(c("navy","white","firebrick3"))(62),
                                      cluster_rows = F,cores = 5,use_gene_short_name=T,show_rownames=T,return_heatmap = T)
ggsave("SF2h.pdf",p,width =6,height = 6)




