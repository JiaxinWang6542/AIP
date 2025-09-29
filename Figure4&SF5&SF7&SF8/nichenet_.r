library(nichenetr)
library(Seurat)
library(tidyverse)
library(circlize)
library(dplyr)
library(clusterProfiler)
library(fgsea)
set.seed(42)
options(stringsAsFactors = F)
setwd("/work/wjx/scRNA-seq/AIP/figure/figure4")


ligand_target_matrix <- readr::read_rds("/work/wjx/scRNA-seq/AIP/codes/cell_communication/NicheNet/database/ligand_target_matrix.rds")
lr_network <- readr::read_rds('/work/wjx/scRNA-seq/AIP/codes/cell_communication/NicheNet/database/lr_network.rds')
weighted_networks <- readr::read_rds('/work/wjx/scRNA-seq/AIP/codes/cell_communication/NicheNet/database/weighted_networks.rds')

lr_network  <- lr_network %>% filter(database !="ppi_prediction_go" & database != "ppi_prediction")

WholeTissueData_B <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")
WholeTissueData_CD4T <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_CD4T_ident.rds.gz")
WholeTissueData_B$celltype2 <- WholeTissueData_B$celltype1
WholeTissueData_CD4T_B <- merge(WholeTissueData_B, WholeTissueData_CD4T)

table(WholeTissueData_CD4T$celltype2)
table(WholeTissueData_CD4T_B$celltype2) 
Idents(WholeTissueData_CD4T_B) <- "celltype2" 





nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = WholeTissueData_CD4T_B, 
                                                receiver_affected = "ABC(IgD-)",
                                                receiver_reference = c("ABC(IgD+)","GC B","Memory B","Naive B","Plasma A","Plasma G"),
                                                sender = c("Tfh","ISG+ T","LIMS1+ Tm","Tem","Temra","TIMP1+ Tm","Tn","TNFRSF18- Treg","TNFRSF18+ Treg"),
                                                ligand_target_matrix = ligand_target_matrix, 
                                                lr_network = lr_network, 
                                                weighted_networks = weighted_networks,
                                                expression_pct = 0.1,
                                                lfc_cutoff = 0.25,
                                                geneset="up",
                                                filter_top_ligands = TRUE, 
                                                top_n_ligands = 10,
                                                top_n_targets = 250, 
                                                cutoff_visualization = 0.33,
                                                verbose = TRUE)



view(nichenet_output$ligand_activities)
nichenet_output$top_ligands
nichenet_output$top_targets
nichenet_output$top_receptors
nichenet_output$ligand_target_matrix
nichenet_output$ligand_target_heatmap
nichenet_output$ligand_target_df
nichenet_output$ligand_expression_dotplot
nichenet_output$ligand_activity_target_heatmap
nichenet_output$ligand_receptor_matrix
nichenet_output$ligand_receptor_heatmap
nichenet_output$ligand_receptor_df
nichenet_output$geneset_oi
nichenet_output$background_expressed_genes


Ligand <- nichenet_output$top_ligands
###ligand activity
ligandActicityRange <- c(min(nichenet_output$ligand_activities[nichenet_output$ligand_activities$test_ligand %in% Ligand,]$aupr_corrected)-0.005,
                         +                          max(nichenet_output$ligand_activities[nichenet_output$ligand_activities$test_ligand %in% Ligand,]$aupr_corrected))
LigandActicityPlot <- ggplot(nichenet_output$ligand_activities,aes(x=1,y=test_ligand))+
  geom_tile(aes(fill=aupr_corrected),color="white",size=0.5)+
  scale_fill_gradient(limit=c(ligandActicityRange),low = "whitesmoke",high = "#ff8c69",breaks=c(0.02,0.025,0.030,0.035))+
  coord_fixed(ratio = 0.8)+
  scale_y_discrete(limit=rev(Ligand))+
  labs(title="Ligand activity",y="Prioritized Ligands",fill="AUPR(target gene prediction ability)")+
  theme(axis.text.x =element_blank(),axis.text.y =element_text(size=13),axis.title.y = element_text(size=10),axis.ticks = element_blank(),axis.title.x = element_blank(),plot.title=element_text(size = 12),legend.position = "bottom")


PlotGenes <- Ligand
Sub_Data <- WholeTissueData_CD4T

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
.cluster_cols <- c(
  "#FB8072", "#1965B0", "#7BAFDE",  "#33A02C", "#B2DF8A","#FF7F00", "#FDB462","#DC050C","#882E72")
Sub_barcolor <- .cluster_cols
names(Sub_barcolor) <- unique(Sub_Data$celltype2)
Sub_SampleBar <-  unique(Sub_Data$celltype2)
names(Sub_SampleBar) <-  unique(Sub_Data$celltype2)


library(ComplexHeatmap)
Sub_ha_bar <- ComplexHeatmap::HeatmapAnnotation(bar = Sub_SampleBar,
                                                col = list(bar = Sub_barcolor),show_legend = F,show_annotation_name = F)

pdf("Figure4h_1.pdf",width = 5.5,height = 5)
ComplexHeatmap::Heatmap(Sub_TF_ClusterMean,cluster_rows = FALSE,cluster_columns = FALSE,cluster_column_slices=F,use_raster = TRUE,name = "Relative expression",
                        col=rev(colorRampPalette(RColorBrewer::brewer.pal(10,"RdYlBu")[1:9],space="rgb")(30)),
                        show_row_names = T,show_column_names = T  ,column_names_rot = 90,top_annotation = Sub_ha_bar)

dev.off()



###ligand and receptor Heatmap
ligand_receptorHeatmap <- nichenet_output$ligand_receptor_matrix[,rev(Ligand)] %>% t %>%
  make_heatmap_ggplot("","Receptors",color = "#B2182B",x_axis_position = "top",legend_title = "Prior interaction potential")+
  theme(axis.text =element_text(size=11))

geneset_oi <- nichenet_output$geneset_oi
active_ligand_target_links_df <- Ligand %>% lapply(get_weighted_ligand_target_links,geneset=geneset_oi,
                                                   ligand_target_matrix=ligand_target_matrix,n=250) %>% bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                  ligand_target_matrix = ligand_target_matrix,cutoff = 0.33)



Order_ligands <- intersect(Ligand,colnames(active_ligand_target_links)) %>% rev() %>% make.names()
Order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()

rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names() 
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()

vis_ligand_target = active_ligand_target_links[Order_targets,Order_ligands] %>% t()


p_ligand_receptor_network <- vis_ligand_target %>% make_heatmap_ggplot("","Predicted target genes",
                                                                       color = "purple",legend_position = "top",
                                                                       x_axis_position = "top",legend_title = "Regulatory potential") +
  theme(axis.text = element_text(face = "italic",size=11)) + scale_fill_gradient2(low = "whitesmoke",  high = "#6a5acd")
p_ligand_receptor_network


targets_defined <- colnames(vis_ligand_target)
targets_defined 
vis_ligand_target_plot <- vis_ligand_target
dim(vis_ligand_target_plot)
p_ligand_receptor_network <- vis_ligand_target_plot %>% make_heatmap_ggplot("","Predicted target genes",
                                                                       color = "purple",legend_position = "top",
                                                                       x_axis_position = "top",legend_title = "Regulatory potential") +
  theme(axis.text = element_text(face = "italic",size=11)) + scale_fill_gradient(low = "#F5F5F5",high = "#9a5acd",breaks = c(0,0.005,0.01))
p_ligand_receptor_network


figures_without_legend <- cowplot::plot_grid(plotlist = list(LigandActicityPlot+theme(legend.position = "none"),
                                                             ligand_receptorHeatmap+theme(legend.position = "none"),
                                                             p_ligand_receptor_network+theme(legend.position = "none")),
                                             ncol = 3,rel_widths = c(2,5.5,10))
Legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(LigandActicityPlot)),
  ggpubr::as_ggplot(ggpubr::get_legend(ligand_receptorHeatmap)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor_network)),
  nrow = 1,
  align = "h",
  rel_widths = c(3,2,2)
)

combined_plot <- cowplot::plot_grid(figures_without_legend,Legends,rel_heights = c(10,5),nrow = 2,align = "hv")
combined_plot
pdf("Figure4h_2.pdf",width = 15,height = 7)
combined_plot
dev.off()



