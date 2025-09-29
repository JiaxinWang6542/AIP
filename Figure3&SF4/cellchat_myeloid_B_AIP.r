setwd("/work/wjx/scRNA-seq/AIP/figure/figure3")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(svglite)
options(stringsAsFactors = FALSE)
WholeTissueData_B=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure2/WholeTissueData_B_ident.rds.gz")
DimPlot(WholeTissueData_B)
WholeTissueData_Myeloid=readRDS("/work/wjx/scRNA-seq/AIP/figure/figure3/WholeTissueData_macrophage_ident.rds.gz")
DimPlot(WholeTissueData_Myeloid)
WholeTissueData_Myeloid_B <- merge(x=WholeTissueData_Myeloid,y=WholeTissueData_B)


WholeTissueData_Myeloid_B_AIP <- subset(WholeTissueData_Myeloid_B,Status == "AIP")
table(WholeTissueData_Myeloid_B_AIP$celltype1)
Idents(WholeTissueData_Myeloid_B_AIP) <- "celltype1"
table(WholeTissueData_Myeloid_B_AIP$Status)
table(WholeTissueData_Myeloid_B_AIP$celltype1)
table(WholeTissueData_Myeloid_B_AIP@active.ident)

data.input  <- WholeTissueData_Myeloid_B_AIP@assays$RNA@data
identity = data.frame(group =WholeTissueData_Myeloid_B_AIP$celltype1, row.names = names(WholeTissueData_Myeloid_B_AIP$celltype1))
unique(identity$group) 

cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "group") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
showDatabaseCategory(CellChatDB)
interaction<- CellChatDB[[1]]

#remove the wrong "CXCL13-CXCR3" ligand-receptor in database
grep("CXCL13_CXCR3",interaction$interaction_name)
CellChatDB[[1]] <- interaction[-710,]
cellchat@DB <- CellChatDB  
unique(CellChatDB$interaction$annotation)
cellchat <- subsetData(cellchat) 
library(future)
plan(multisession,workers=20)
options(future.globals.maxSize= 2097152000)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.1,raw.use = T,population.size = F)

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)


vertex.receiver = seq(1,4)
netVisual_heatmap(cellchat, color.heatmap = "Reds",signaling = "CCL",sources.use = unique(WholeTissueData_Myeloid$celltype1),remove.isolate = T)
pdf("Figure3f.pdf",width = 6,height = 6)
netVisual_aggregate(cellchat,  signaling.name = "CXCL/CCL",layout = "circle",signaling = c("CCL","CXCL") ,sources.use = unique(WholeTissueData_Myeloid$celltype1),targets.use =unique(WholeTissueData_B$celltype1),top=0.2,edge.width.max = 16)
dev.off()

cellchat@netP$pathways
pathways.show <- c("CCL","CXCL") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle")
pdf("Figure3g.pdf",width = 8,height = 5)
netVisual_chord_gene(cellchat, sources.use = unique(WholeTissueData_Myeloid$celltype1), targets.use = "ABC(IgD-)", signaling = c("CCL","CXCL"),legend.pos.x = 8)
dev.off()

