setwd("/work/wjx/scRNA-seq/AIP")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(svglite)
options(stringsAsFactors = FALSE)
WholeTissueData_B=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")
DimPlot(WholeTissueData_B)
WholeTissueData_CD4T=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_CD4T_ident.rds.gz")
DimPlot(WholeTissueData_CD4T)
WholeTissueData_CD4T$celltype1 <- WholeTissueData_CD4T$celltype2
WholeTissueData_CD4T_B <- merge(x=WholeTissueData_CD4T,y=WholeTissueData_B)
table(WholeTissueData_CD4T_B$celltype1)


WholeTissueData_CD4T_B_AIP <- subset(WholeTissueData_CD4T_B,Status == "AIP")
table(WholeTissueData_CD4T_B_AIP$celltype1)
Idents(WholeTissueData_CD4T_B_AIP) <- "celltype1"
table(WholeTissueData_CD4T_B_AIP$Status)
table(WholeTissueData_CD4T_B_AIP$celltype1)
table(WholeTissueData_CD4T_B_AIP@active.ident)


data.input  <- WholeTissueData_CD4T_B_AIP@assays$RNA@data
identity = data.frame(group =WholeTissueData_CD4T_B_AIP$celltype1, row.names = names(WholeTissueData_CD4T_B_AIP$celltype1)) # create a dataframe consisting of the cell labels
unique(identity$group) 

cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "group") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
showDatabaseCategory(CellChatDB)

cellchat@DB <- CellChatDB  
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat)
library(future)
plan(multisession,workers=20)
options(future.globals.maxSize= 2097152000)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.07,raw.use = T,population.size = F)

cellchat <- filterCommunication(cellchat, min.cells = 10)



cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)

saveRDS(cellchat,file = "Cellchat_CD4T_B_AIP_20240126",compress = "gzip")


#######Normal cellchat
WholeTissueData_CD4T_B_Normal <- subset(WholeTissueData_CD4T_B,Status == "Normal")
table(WholeTissueData_CD4T_B_Normal$celltype1)
Idents(WholeTissueData_CD4T_B_Normal) <- "celltype1"
table(WholeTissueData_CD4T_B_Normal$Status)
table(WholeTissueData_CD4T_B_Normal$celltype1)
table(WholeTissueData_CD4T_B_Normal@active.ident)


data.input  <- WholeTissueData_CD4T_B_Normal@assays$RNA@data
identity = data.frame(group =WholeTissueData_CD4T_B_Normal$celltype1, row.names = names(WholeTissueData_CD4T_B_Normal$celltype1)) # create a dataframe consisting of the cell labels
unique(identity$group) 

cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "group") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
showDatabaseCategory(CellChatDB)

cellchat@DB <- CellChatDB  
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) 
library(future)
plan(multisession,workers=20)
options(future.globals.maxSize= 2097152000)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.07,raw.use = T,population.size = F)

cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)

saveRDS(cellchat,file = "Cellchat_CD4T_B_Normal_20240126",compress = "gzip")






cellchat.CD4T.B.AIP <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/result/cellchat/Cellchat_CD4T_B_AIP_20240126")
levels(cellchat.CD4T.B.AIP@idents)
cellchat.CD4T.B.N <- readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/result/cellchat/Cellchat_CD4T_B_Normal_20240126")
levels(cellchat.CD4T.B.N@idents)

group.new = levels(cellchat.AIP@idents)
cellchat.N <- liftCellChat(cellchat.N, group.new)
object.list <- list(N = cellchat.N, AIP = cellchat.AIP)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf("Figure4g.pdf",width = 6,height = 5.5)
gg2
dev.off()


a <- netVisual_bubble(cellchat.CD4T.B.AIP, sources.use = "ABC(IgD-)", targets.use = c("Tfh"), angle.x = 45, return.data = T)
a <- a[[1]]
Sub_Data <- top_n(a,n=15,wt=prob)
pdf("SF5h_1.pdf",width = 3.5,height = 7.5)
ggplot(Sub_Data,aes(x=source.target,y=interaction_name_2))+
  geom_point(aes(color=prob),size=5)+
  scale_color_gradientn(colors = colorRampPalette(c(RColorBrewer::brewer.pal(6,'Reds')),space='rgb')(10),limits=c(0.1,0.5))+
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=90,vjust = 0.5,hjust = 1,color="black"),
        axis.text.y=element_text(size=12,color="black"),axis.title=element_blank(),legend.position = "right")
dev.off()

b <- netVisual_bubble(cellchat.CD4T.B.AIP, sources.use = c("Tfh"), targets.use = c("ABC(IgD-)"), angle.x = 45, return.data = T)
b <- b[[1]]
Sub_Data <- top_n(b,n=15,wt=prob)
pdf("SF5h_2.pdf",width = 3.5,height = 7.5)
ggplot(Sub_Data,aes(x=source.target,y=interaction_name_2))+
  geom_point(aes(color=prob),size=5)+
  scale_color_gradientn(colors = colorRampPalette(c(RColorBrewer::brewer.pal(6,'Reds')),space='rgb')(10),limits=c(0.1,0.5))+
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=90,vjust = 0.5,hjust = 1,color="black"),
        axis.text.y=element_text(size=12,color="black"),axis.title=element_blank(),legend.position = "right")
dev.off()






######CD8T
WholeTissueData_B=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")
DimPlot(WholeTissueData_B)
WholeTissueData_CD8T=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure4/WholeTissueData_CD8T_ident.rds.gz")
DimPlot(WholeTissueData_CD8T)
WholeTissueData_CD8T$celltype1 <- WholeTissueData_CD8T$celltype2
table(WholeTissueData_CD8T$celltype1)
WholeTissueData_CD8T_B <- merge(x=WholeTissueData_CD8T,y=WholeTissueData_B)
table(WholeTissueData_CD8T_B$celltype1)


WholeTissueData_CD8T_B_AIP <- subset(WholeTissueData_CD8T_B,Status == "AIP")
table(WholeTissueData_CD8T_B_AIP$celltype1)
Idents(WholeTissueData_CD8T_B_AIP) <- "celltype1"
table(WholeTissueData_CD8T_B_AIP$Status)
table(WholeTissueData_CD8T_B_AIP$celltype1)
table(WholeTissueData_CD8T_B_AIP@active.ident)


data.input  <- WholeTissueData_CD8T_B_AIP@assays$RNA@data
identity = data.frame(group =WholeTissueData_CD8T_B_AIP$celltype1, row.names = names(WholeTissueData_CD8T_B_AIP$celltype1)) # create a dataframe consisting of the cell labels
unique(identity$group) 

cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "group",levels = c(levels(WholeTissueData_B$celltype1),"Trm","Temra","MAIT","Tem","Tn")) # set "group" as default cell identity
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
showDatabaseCategory(CellChatDB)

cellchat@DB <- CellChatDB  
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) 
library(future)
plan(multisession,workers=20)
options(future.globals.maxSize= 2097152000)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.07,raw.use = T,population.size = F)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)


pdf("SF7e.pdf",width = 5,height = 4)
netVisual_heatmap(cellchat, measure="weight",sources.use=c("Trm","Temra","MAIT","Tem","Tn"),targets.use =c("ABC(IgD-)"),color.heatmap = "YlOrRd",remove.isolate = T)
dev.off()



a <- netVisual_bubble(cellchat, sources.use = "Trm", targets.use = "ABC(IgD-)", angle.x = 45, return.data = T)
a <- a[[1]]
Sub_Data <- top_n(a,n=15,wt=prob) %>% arrange(-prob)
pdf("SF7f.pdf",width = 7.5,height = 3)
ggplot(Sub_Data,aes(x=source.target,y=interaction_name_2))+
  geom_point(aes(color=prob),size=5)+
  scale_color_gradientn(colors = colorRampPalette(c(RColorBrewer::brewer.pal(6,'Reds')),space='rgb')(10),limits=c(0.1,0.5))+
  scale_y_discrete(limits = unique(Sub_Data$interaction_name_2))+
  coord_flip()+
  theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),axis.text.x=element_text(size=12,angle=45,vjust = 1,hjust = 1,color="black"),
        axis.text.y=element_text(size=12,color="black"),axis.title=element_blank(),legend.position = "right")
dev.off()


