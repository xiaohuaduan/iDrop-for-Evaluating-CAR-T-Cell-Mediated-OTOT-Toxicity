suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(patchwork)
})


ptm = Sys.time()
Allo_Control<- subset(Colon_T, subset = sample_group %in% "Allo_Control")
Allo_Utr_T<- subset(Colon_T, subset = sample_group %in% "Allo_Utr_T")
Allo_Her2_T<- subset(Colon_T, subset = sample_group %in% "Allo_Her2_T")
Auto_Control<- subset(Colon_T, subset = sample_group %in% "Auto_Control")
Auto_Utr_T<- subset(Colon_T, subset = sample_group %in% "Auto_Utr_T")
Auto_Her2_T<- subset(Colon_T, subset = sample_group %in% "Auto_Her2_T")

Allo_Control_chat <- createCellChat(object = Allo_Control, group.by = "clusters", assay = "RNA")
Allo_Utr_T_chat <- createCellChat(object = Allo_Utr_T, group.by = "clusters", assay = "RNA")
Allo_Her2_T_chat <- createCellChat(object = Allo_Her2_T, group.by = "clusters", assay = "RNA")
Auto_Control_chat <- createCellChat(object = Auto_Control, group.by = "clusters", assay = "RNA")
Auto_Utr_T_chat <- createCellChat(object = Auto_Utr_T, group.by = "clusters", assay = "RNA")
Auto_Her2_T_chat <- createCellChat(object = Auto_Her2_T, group.by = "clusters", assay = "RNA")

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)
Allo_Control_chat@DB <- CellChatDB.use
Allo_Utr_T_chat@DB <- CellChatDB.use
Allo_Her2_T_chat@DB <- CellChatDB.use
Auto_Control_chat@DB <- CellChatDB.use
Auto_Utr_T_chat@DB <- CellChatDB.use
Auto_Her2_T_chat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
Allo_Control_chat <- subsetData(Allo_Control_chat) # This step is necessary even if using the whole database

Allo_Control_chat <- identifyOverExpressedGenes(Allo_Control_chat)
Allo_Control_chat <- identifyOverExpressedInteractions(Allo_Control_chat)

Allo_Utr_T_chat <- subsetData(Allo_Utr_T_chat) # This step is necessary even if using the whole database

Allo_Utr_T_chat <- identifyOverExpressedGenes(Allo_Utr_T_chat)
Allo_Utr_T_chat <- identifyOverExpressedInteractions(Allo_Utr_T_chat)

Allo_Her2_T_chat <- subsetData(Allo_Her2_T_chat) # This step is necessary even if using the whole database

Allo_Her2_T_chat <- identifyOverExpressedGenes(Allo_Her2_T_chat)
Allo_Her2_T_chat <- identifyOverExpressedInteractions(Allo_Her2_T_chat)

Auto_Control_chat <- subsetData(Auto_Control_chat) # This step is necessary even if using the whole database

Auto_Control_chat <- identifyOverExpressedGenes(Auto_Control_chat)
Auto_Control_chat <- identifyOverExpressedInteractions(Auto_Control_chat)

Auto_Utr_T_chat <- subsetData(Auto_Utr_T_chat) # This step is necessary even if using the whole database

Auto_Utr_T_chat <- identifyOverExpressedGenes(Auto_Utr_T_chat)
Auto_Utr_T_chat <- identifyOverExpressedInteractions(Auto_Utr_T_chat)

Auto_Her2_T_chat <- subsetData(Auto_Her2_T_chat) # This step is necessary even if using the whole database

Auto_Her2_T_chat <- identifyOverExpressedGenes(Auto_Her2_T_chat)
Auto_Her2_T_chat <- identifyOverExpressedInteractions(Auto_Her2_T_chat)

#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)

Allo_Control_chat  <- smoothData(Allo_Control_chat , adj = PPI.human)
Allo_Utr_T_chat  <- smoothData(Allo_Utr_T_chat , adj = PPI.human)
Allo_Her2_T_chat  <- smoothData(Allo_Her2_T_chat , adj = PPI.human)
Auto_Control_chat  <- smoothData(Auto_Control_chat , adj = PPI.human)
Auto_Utr_T_chat  <- smoothData(Auto_Utr_T_chat , adj = PPI.human)
Auto_Her2_T_chat  <- smoothData(Auto_Her2_T_chat , adj = PPI.human)

ptm = Sys.time()
Allo_Control_chat <- computeCommunProb(Allo_Control_chat, type =  "truncatedMean", trim = 0.1,population.size = TRUE,nboot = 20)
Allo_Utr_T_chat <- computeCommunProb(Allo_Utr_T_chat, type =  "truncatedMean", trim = 0.1,population.size = TRUE,nboot = 20)
Allo_Her2_T_chat <- computeCommunProb(Allo_Her2_T_chat, type =  "truncatedMean", trim = 0.1,population.size = TRUE,nboot = 20)
Auto_Control_chat <- computeCommunProb(Auto_Control_chat, type =  "truncatedMean", trim = 0.1,population.size = TRUE,nboot = 20)
Auto_Utr_T_chat <- computeCommunProb(Auto_Utr_T_chat, type =  "truncatedMean", trim = 0.1,population.size = TRUE,nboot = 20)
Auto_Her2_T_chat <- computeCommunProb(Auto_Her2_T_chat, type =  "truncatedMean", trim = 0.1,population.size = TRUE,nboot = 20)

Allo_Control_chat <- filterCommunication(Allo_Control_chat, min.cells = 10)
Allo_Utr_T_chat <- filterCommunication(Allo_Utr_T_chat, min.cells = 10)
Allo_Her2_T_chat <- filterCommunication(Allo_Her2_T_chat, min.cells = 10)
Auto_Control_chat <- filterCommunication(Auto_Control_chat, min.cells = 10)
Auto_Utr_T_chat <- filterCommunication(Auto_Utr_T_chat, min.cells = 10)
Auto_Her2_T_chat <- filterCommunication(Auto_Her2_T_chat, min.cells = 10)

Allo_Control_chat <- computeCommunProbPathway(Allo_Control_chat)
Allo_Utr_T_chat <- computeCommunProbPathway(Allo_Utr_T_chat)
Allo_Her2_T_chat <- computeCommunProbPathway(Allo_Her2_T_chat)
Auto_Control_chat <- computeCommunProbPathway(Auto_Control_chat)
Auto_Utr_T_chat <- computeCommunProbPathway(Auto_Utr_T_chat)
Auto_Her2_T_chat <- computeCommunProbPathway(Auto_Her2_T_chat)

Allo_Control_chat <- aggregateNet(Allo_Control_chat)
Allo_Utr_T_chat <- aggregateNet(Allo_Utr_T_chat)
Allo_Her2_T_chat <- aggregateNet(Allo_Her2_T_chat)
Auto_Control_chat <- aggregateNet(Auto_Control_chat)
Auto_Utr_T_chat <- aggregateNet(Auto_Utr_T_chat)
Auto_Her2_T_chat <- aggregateNet(Auto_Her2_T_chat)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
par(mfrow = c(2,3), xpd=TRUE)
groupSize <- as.numeric(table(Allo_Control_chat@idents))

netVisual_circle(Allo_Control_chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Allo_Control_chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ptm = Sys.time()
groupSize <- as.numeric(table(Allo_Utr_T_chat@idents))

netVisual_circle(Allo_Utr_T_chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Allo_Utr_T_chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

ptm = Sys.time()
groupSize <- as.numeric(table(Allo_Her2_T_chat@idents))

netVisual_circle(Allo_Her2_T_chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Allo_Her2_T_chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

groupSize <- as.numeric(table(Auto_Control_chat@idents))

netVisual_circle(Auto_Control_chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Auto_Control_chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

groupSize <- as.numeric(table(Auto_Utr_T_chat@idents))

netVisual_circle(Auto_Utr_T_chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Auto_Utr_T_chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

groupSize <- as.numeric(table(Auto_Her2_T_chat@idents))

netVisual_circle(Auto_Her2_T_chat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(Auto_Her2_T_chat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



## --- Per-sender circle plots with correct vertex weights/order ---

# 1) Start from your matrix
M <- Allo_Her2_T_chat@net$weight
stopifnot(is.matrix(M))

# 2) Keep only groups present as BOTH senders and receivers, same order for rows/cols
common <- intersect(rownames(M), colnames(M))
if (length(common) < 2) stop("Too few common groups in the matrix.")
M <- M[common, common, drop = FALSE]

# 3) Derive vertex weights (group sizes) in the SAME order as matrix rows
#    Use the CellChat idents to count cells per group
gs_tbl <- table(factor(Allo_Control_chat@idents, levels = common))
groupSize <- as.numeric(gs_tbl)           # length == nrow(M)
names(groupSize) <- common

# (Optional) if all group sizes are zero (shouldn’t happen), fall back to row sums
if (all(groupSize == 0)) {
  groupSize <- as.numeric(rowSums(M, na.rm = TRUE))
  names(groupSize) <- common
}

# 4) Global max for edge scaling
wmax <- max(M, na.rm = TRUE)

# 5) Plot: outgoing edges from each sender i
op <- par(mfrow = c(2, 4), xpd = TRUE)   # your layout
on.exit(par(op), add = TRUE)

for (i in seq_len(nrow(M))) {
  mat2 <- matrix(0, nrow = nrow(M), ncol = ncol(M), dimnames = dimnames(M))
  v <- M[i, ]
  v[!is.finite(v)] <- 0
  mat2[i, ] <- v
  
  netVisual_circle(
    mat2,
    vertex.weight   = groupSize[rownames(mat2)],  # exact vertex order
    weight.scale    = TRUE,
    label.edge      = FALSE,
    vertex.label.cex= 0.9,
    edge.weight.max = wmax,
    title.name      = paste0(rownames(M)[i], " → all (weight)")
  )
}





par(mfrow = c(2,2), xpd=TRUE)
netVisual_bubble(Allo_Utr_T_chat, sources.use = c(1), targets.use = c(6), remove.isolate = T)
netVisual_bubble(Allo_Her2_T_chat, sources.use =  c(1), targets.use = c(6), remove.isolate = T)
netVisual_bubble(Auto_Utr_T_chat, sources.use =  c(1), targets.use = c(6), remove.isolate = T)
netVisual_bubble(Auto_Her2_T_chat, sources.use =  c(1), targets.use = c(6), remove.isolate = T)


object.list <- list(NL1 = Allo_Utr_T_chat, LS1 = Allo_Her2_T_chat, NL2 = Auto_Utr_T_chat, LS2 = Auto_Her2_T_chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_bubble(cellchat, sources.use =  c(1), targets.use = c(6), comparison = c(3,4,1,2),remove.isolate = T, max.dataset = c(1,2))




