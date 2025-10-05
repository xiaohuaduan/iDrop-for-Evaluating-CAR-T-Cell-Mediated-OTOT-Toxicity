library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(harmony)
library(ggrepel)
library(dplyr)
# Creat seurat objects
sample_paths <- list(
  Auto_Control = "D:/new data/xilis/combine singel cell dataset/control/raw_matrix",
  Auto_Utr_T = "D:/new data/xilis/combine singel cell dataset/utr/raw_matrix",
  Auto_Her2_T = "D:/new data/xilis/combine singel cell dataset/her2/raw_matrix",
  Allo_Control = "D:/new data/xilis/combine singel cell dataset/raw_matrix_h1_colon_control",
  Allo_Utr_T = "D:/new data/xilis/combine singel cell dataset/colon_utr_raw_matrix",
  Allo_Her2_T = "D:/new data/xilis/combine singel cell dataset/colon_her2_raw_matrix"
)
seurat_list <- list()
# Load, QC, Normalize, DoubletFinder
for (sample_name in names(sample_paths)) {
  cat("Processing:", sample_name, "\n")
  
  # Read 10X
  counts <- Read10X(data.dir = sample_paths[[sample_name]])
  obj <- CreateSeuratObject(counts = counts, project = sample_name)
  obj$orig.ident <- sample_name
  
  # QC
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 250 & nFeature_RNA < 7500 & percent.mt < 5)
  
  # Normalize, variable genes, PCA
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(obj, PCs = 1:20,sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  ggplot(bcmvn,aes(pK,BCmetric,group=1))+geom_point()+geom_line()
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  ## Homotypic doublet proportion estimate
  annotations <- obj@meta.data$orig.ident
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(obj@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  # run DoubletFinder
  obj<- doubletFinder(seu = obj, 
                      PCs = 1:20, 
                      pK = optimal.pk,
                      nExp = nExp.poi.adj)
  metadata <- obj@meta.data
  colnames(metadata)[6] <- "doublet_finder"
  obj@meta.data <- metadata 
  
  # subset and save
  obj <- subset(obj, doublet_finder == "Singlet")
  
  seurat_list[[sample_name]] <- obj
  
}

# Merge all samples into one Seurat object harmony
Colon_T<- merge(x =seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "Merged")
Colon_T$sample_group <- Idents(Colon_T)
# Normalize, variable genes, scale, PCA again on merged object
Colon_T <- NormalizeData(Colon_T)
Colon_T <- FindVariableFeatures(Colon_T, selection.method = "vst", nfeatures = 3000)
Colon_T <- ScaleData(Colon_T, verbose = FALSE)

Colon_T <- RunPCA(Colon_T, npcs = 50, verbose = FALSE)
ElbowPlot(Colon_T, ndims = ncol(Embeddings(Colon_T, "pca")))
PCHeatmap(Colon_T, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)

# Run Harmony
Colon_T <- RunHarmony(Colon_T, group.by.vars = "orig.ident",dims.use = 1:30, max.iter.harmony = 50)

# UMAP, clustering using Harmony-corrected PCA
Colon_T <- RunUMAP(Colon_T, reduction = "harmony", dims = 1:30)
Colon_T <- FindNeighbors(Colon_T, reduction = "harmony", dims = 1:30)
Colon_T <- FindClusters(Colon_T, resolution = 1)


# Visualization
DimPlot(Colon_T, reduction = "umap", label = TRUE)

DimPlot(Colon_T, reduction = "umap", label = TRUE,split.by = "sample_group")


# Check
table(Colon_T$sample_group)
table(Colon_T$orig.ident)

# FindAllMarkers, heatmap, save csv.file
Colon_T<- JoinLayers(Colon_T)
Colon_T_markers <- FindAllMarkers(Colon_T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
Colon_T_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10_Colon_T_markers <- Colon_T_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Colon_T, features = top10_Colon_T_markers$gene) + NoLegend()
write.csv(Colon_T_markers, "D:/new data/xilis/combine singel cell dataset/Colon_T_markers.csv",row.names = T)

# save dataset
saveRDS(Colon_T, file="Colon_T.rds")
Colon_T<-readRDS("Colon_T.rds")

#  ---------- Visualization ----------


Idents(Colon_T) <- "orig.ident"
Colon_T$Samples <- Idents(Colon_T)

sample_colors <- c(
  "Allo_Control" = "#0a9396",   # teal
  "Allo_Utr_T"   = "#ee9b00",   # orange
  "Allo_Her2_T"  = "#ae2012",   # red
  "Auto_Control" = "#005f73",   # dark teal
  "Auto_Utr_T"   = "#FF8C00",   # dark orange
  "Auto_Her2_T"  = "#8B0000"   # dark red
)
DimPlot(
  Colon_T,
  group.by = "Samples",
  label = F,
  cols = sample_colors,
  pt.size = 0.5
)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 20))

DimPlot(
  Colon_T,
  split.by = "Samples",
  label = F,
  cols = sample_colors,
  pt.size = 0.5
)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 20))

# cell anotation 

Colon_T <- FindNeighbors(Colon_T, reduction = "harmony", dims = 1:30)
Colon_T <- FindClusters(Colon_T, resolution = 1.5)
Colon_T <- RunUMAP(Colon_T, reduction = "harmony",dims = 1:30)
Colon_T_ident <- setNames(c("Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop",
                            "Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop",
                            "Colon_iDrop","Colon_iDrop","Colon_iDrop","T_cell","Colon_iDrop","Colon_iDrop","Colon_iDrop",
                            "Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","T_cell","Colon_iDrop","T_cell",
                            "Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop"),
                          levels(Colon_T))
Colon_T<- RenameIdents(Colon_T, Colon_T_ident)
Colon_T$Immu_iDrop <- Idents(Colon_T)
DimPlot(Colon_T, reduction = "umap", label = T)+NoLegend()+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 20))

DimPlot(Colon_T, reduction = "umap", label = TRUE,split.by = "orig.ident",repel = T)

Colon_T <- FindClusters(Colon_T, resolution = 1.5)
Colon_T <- RunUMAP(Colon_T, reduction = "harmony",dims = 1:30)
DimPlot(Colon_T, reduction = "umap", label = T)+NoLegend()+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 20))
DotPlot(Colon_T,c("MAP2","LGR5","FABP1","KRT20","EPHB2","MUC5AC","ERBB2","COL1A1","TOP2A"))

Colon_T_ident1 <- setNames(c("Enterocyte","Enterocyte","Goblet cell","Fibroblast","Enterocyte","Neuronal cell","Goblet cell",
                            "Enterocyte","Enterocyte","Enterocyte","Fibroblast","Unknown","Neuronal cell","Goblet cell",
                            "Fibroblast","T cell","Neuronal cell","Neuronal cell","Neuronal cell","Neuronal cell",
                            "Enteroendocrine cell","Neuronal cell","Neuronal cell","T cell",
                            "Fibroblast","T cell","Enteroendocrine cell","Enteroendocrine cell","Poliferating cell","Enterocyte","Neuronal cell","Goblet cell"),
                            levels(Colon_T))

Colon_T<- RenameIdents(Colon_T, Colon_T_ident1)
Colon_T$clusters <- Idents(Colon_T)
DimPlot(Colon_T, reduction = "umap", label = T)+NoLegend()

DimPlot(Colon_T, reduction = "umap", label = TRUE,split.by = "orig.ident",repel = T)


#T cells
library(dplyr)
Auto <- subset(Colon_T, subset = sample_group %in% c("Auto_Control", "Auto_Utr_T","Auto_Her2_T"))

DimPlot(Auto, reduction = "umap", label = T)+NoLegend()+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 20))


DimPlot(Auto, reduction = "umap", label = TRUE,split.by = "orig.ident",repel = T)


Idents(Colon_T) <- "orig.ident"
Colon_T$Samples <- Idents(Colon_T)

T_cell<- subset(Colon_T, subset = Immu_iDrop == "T_cell")

Idents(T_cell) <- "Samples"

DimPlot(T_cell, reduction = "umap", label = F)

de_Her2_vs_Untr_T <- FindMarkers(
  T_cell,
  ident.1 = "Auto_Her2_T",
  ident.2 = "Auto_Utr_T",
  logfc.threshold = 0.25,
  test.use = "wilcox"  # default
)

sig_de_Her2_vs_Untr_T <- de_Her2_vs_Untr_T %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

head(sig_de_Her2_vs_Untr_T)

write.csv(sig_de_Her2_vs_Untr_T, "D:/new data/xilis/combine singel cell dataset/sig_de_Her2_vs_Untr_T.csv",row.names = T)


de_Her2_vs_Untr_T$gene <- rownames(de_Her2_vs_Untr_T)

de_Her2_vs_Untr_T$sig <- with(de_Her2_vs_Untr_T,
                              ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, "Significant", "Not"))

# Example: adjust this to your DE result
res_df <- de_Her2_vs_Untr_T
res_df$gene <- rownames(res_df)

# Label significant genes
res_df <- res_df %>%
  mutate(
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC > 0.5 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Count numbers
n_up <- sum(res_df$sig == "Up")
n_down <- sum(res_df$sig == "Down")

# Top genes to label
top_genes <- res_df %>%
  filter(sig != "NS") %>%
  arrange(p_val_adj) %>%
  slice(1:15)

# Plot
ggplot(res_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = sig), alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "T cell: Her2_T vs Utr_T",
    x = "logFC",
    y = "-log10(p)"
  ) +
  theme(legend.position = "top")

# T cell activation pathway

sample_order <- c(
  "Auto_Utr_T",
  "Allo_Utr_T",
  "Auto_Her2_T",
  "Allo_Her2_T"
)
T_cell$Samples <- factor(T_cell$Samples, levels = sample_order)

DotPlot(T_cell, c("TNFSF10","GZMB","GZMA","CSF1","CSF2","IL13","IL5","EIF2AK2",
                     "ICAM1","OAS1","CXCL10","CXCL9","CXCL11","CXCL5","XCL1"),group.by = "Samples")+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y  = element_text(face = "italic"))+coord_flip()

table(T_cell$Samples)


## ---------- Setup ----------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

obj <- T_cell                           # your Seurat object
assay_use <- DefaultAssay(obj)          # e.g. "RNA"

## Enforce your sample order (adjust the column name if needed)
sample_order <- c("Auto_Utr_T","Allo_Utr_T","Auto_Her2_T","Allo_Her2_T")
obj$Samples <- factor(obj$Samples, levels = sample_order)

## ---------- 1) Define gene sets ----------
gene_sets <- list(
  Extrinsic_and_intrinsic_apoptosis_inducers = c("GZMA","GZMB","TNFSF10","LTA","TNF"),
  Cytokine_and_inflammatory_response         = c("IL5","IL13","CSF1","CSF2"),
  Type_II_interferon_signaling               = c("EIF2AK2","ICAM1","OAS1","CXCL10","CXCL9","ISG15","IFNG"),
  Chemokine_Signaling                        = c("CXCL9","CXCL10","CXCL11","CXCL5","XCL1")
)

## ---------- 2) Keep only genes present ----------
present <- rownames(obj[[assay_use]])
features_use <- lapply(gene_sets, function(gs) intersect(gs, present))

lens <- vapply(features_use, length, 1L)
if (any(lens < 3)) {
  warning("These sets have <3 genes present and may yield unstable scores: ",
          paste(names(lens)[lens < 3], collapse = ", "))
}

## ---------- 3) AddModuleScore ----------
obj <- AddModuleScore(
  object  = obj,
  features = features_use,
  name     = "SIG_",
  assay    = assay_use,
  nbin = 24, ctrl = 100, seed = 1
)

## ---------- 4) Rename score columns to readable labels ----------
n_sets  <- length(features_use)
old_cols <- paste0("SIG_", seq_len(n_sets))
old_cols <- old_cols[old_cols %in% colnames(obj@meta.data)]
stopifnot(length(old_cols) == n_sets)

new_cols <- names(features_use)
colnames(obj@meta.data)[match(old_cols, colnames(obj@meta.data))] <- new_cols

## ---------- 5) DotPlot-style visualization (uses AddModuleScore outputs) ----------
pathways_to_plot <- c(
  "Extrinsic_and_intrinsic_apoptosis_inducers",
  "Cytokine_and_inflammatory_response",
  "Type_II_interferon_signaling",
  "Chemokine_Signaling"
)
pathways_to_plot <- intersect(pathways_to_plot, colnames(obj@meta.data))
stopifnot(length(pathways_to_plot) == 4)

# summarize: mean score (color) & % cells with score > 0 (size)
dot_df <- obj@meta.data %>%
  dplyr::select(Samples, all_of(pathways_to_plot)) %>%
  tidyr::pivot_longer(all_of(pathways_to_plot), names_to = "Pathway", values_to = "Score") %>%
  dplyr::group_by(Samples, Pathway) %>%
  dplyr::summarise(mean_score = mean(Score, na.rm = TRUE),
                   pct_pos    = mean(Score > 0, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::group_by(Pathway) %>%
  dplyr::mutate(z_mean = scale(mean_score)[,1]) %>%
  dplyr::ungroup()

dot_df$Samples <- factor(dot_df$Samples, levels = sample_order)
dot_df$Pathway <- factor(dot_df$Pathway, levels = pathways_to_plot)

ggplot(dot_df, aes(x = Samples, y = Pathway)) +
  geom_point(aes(size = pct_pos, color = z_mean)) +
  scale_size(range = c(2.5, 9), name = "Pct > 0") +
  scale_color_viridis_c(name = "Mean score (z)", option = "D") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module scores (AddModuleScore) — DotPlot-style by sample", x = NULL, y = NULL)


#colon cells

Colon <- subset(
  Colon_T,
  idents = c("T cell"),
  invert = TRUE
)

Idents(Colon) <- "clusters"
Colon@active.ident <- factor(Colon@active.ident, 
                            levels=c("Enterocyte",
                                     "Goblet cell",
                                     "Unknown",
                                     "Fibroblast",
                                     "Enteroendocrine cell",
                                     "Poliferating cell",
                                     "Neuronal cell"))
DimPlot(Colon, pt.size = 0.5)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 14, face = "bold"))


# Define your color vector
colors16 <- c(
  "#0a9396",  # sea green
  "#005f73",  # dark teal
  "#a6cee3",  # light blue (Set3)
  "#1f78b4",  # blue (Set1)
  "#b2df8a",  # light green
  "#33a02c",  # green
  "#fb9a99",  # pink
  "#e31a1c",  # red
  "#fdbf6f",  # light orange
  "#ff7f00",  # orange
  "#cab2d6",  # lavender
  "#6a3d9a",  # purple
  "#ffff99",  # pale yellow
  "#b15928",  # brown
  "#8e9aaf",   # slate gray (added 16th to balance warm/cool)
  "#001219"  # deep navy
)
# Ensure factor levels for consistent order

Colon$samples <- factor(Colon$orig.ident, levels = c("Allo_Control", "Allo_Utr_T","Allo_Her2_T"))

# Assign to cluster identities
names(colors16) <- levels(Colon)

DimPlot(Colon,cols = colors16,pt.size = 0.5)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 14))

DimPlot(Colon, pt.size = 0.1,split.by = "Samples",ncol = 3,cols = colors16)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 14))

DotPlot(Colon,  c("FABP1","MUC5AC","LUM","CHGA","MKI67","PAX6"))+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y  = element_text(face = "italic"))+coord_flip()
DotPlot(Auto_Colon,  c("FABP1","MUC5AC","COL1A1","PAX6","CHGA","MKI67"))+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y  = element_text(face = "italic"))+coord_flip()

DotPlot(Colon, c("BCL2","RIPK1","BAK1","NFKBIA",
                         "PMAIP1","FAS","FADD",
                         "TNFRSF1B","CASP8",
                         "CASP7","CASP4","CASP1","IRF7",
                         "IRF3","IRF2","IRF1","BIRC3"),group.by = "Samples")+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y = element_text(face = "italic"))+coord_flip()

# colon iDrop apoptosis pathway

sample_order <- c(
  "Auto_Control",
  "Auto_Utr_T",
  "Auto_Her2_T",
  "Allo_Control",
  "Allo_Utr_T",
  "Allo_Her2_T"
)
Colon$Samples <- factor(Colon$Samples, levels = sample_order)

DotPlot(Colon, c("BCL2","RIPK1","BAK1","BBC3","NFKBIA",
                 "PMAIP1","FAS","FADD","TNFRSF10B",
                 "TNFRSF1B","CASP8",
                 "CASP7","CASP4","CASP1","IRF7",
                 "IRF3","IRF2","IRF1","BIRC3"),group.by = "Samples")+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y = element_text(face = "italic"))+coord_flip()

table(T_cell$Samples)


## ---------- Setup ----------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

obj <- Colon                           # your Seurat object
assay_use <- DefaultAssay(obj)          # e.g. "RNA"

## Enforce your sample order (adjust the column name if needed)
sample_order <- c(
  "Auto_Control",
  "Auto_Utr_T",
  "Auto_Her2_T",
  "Allo_Control",
  "Allo_Utr_T",
  "Allo_Her2_T"
)
obj$Samples <- factor(obj$Samples, levels = sample_order)

## ---------- 1) Define gene sets ----------
gene_sets <- list(
  apoptosis_pathway = c("BCL2","RIPK1","BAK1","BBC3","NFKBIA",
    "PMAIP1","FAS","FADD","TNFRSF10B",
    "TNFRSF1B","CASP8",
    "CASP7","CASP4","CASP1","IRF7",
    "IRF3","IRF2","IRF1","BIRC3")
)

## ---------- 2) Keep only genes present ----------
present <- rownames(obj[[assay_use]])
features_use <- lapply(gene_sets, function(gs) intersect(gs, present))

lens <- vapply(features_use, length, 1L)
if (any(lens < 3)) {
  warning("These sets have <3 genes present and may yield unstable scores: ",
          paste(names(lens)[lens < 3], collapse = ", "))
}

## ---------- 3) AddModuleScore ----------
obj <- AddModuleScore(object  = obj,
  features = features_use,
  name     = "SIG_",
  assay    = assay_use,
  nbin = 24, ctrl = 100, seed = 1
)

## ---------- 4) Rename score columns to readable labels ----------
n_sets  <- length(features_use)
old_cols <- paste0("SIG_", seq_len(n_sets))
old_cols <- old_cols[old_cols %in% colnames(obj@meta.data)]
stopifnot(length(old_cols) == n_sets)

new_cols <- names(features_use)
colnames(obj@meta.data)[match(old_cols, colnames(obj@meta.data))] <- new_cols

## ---------- 5) DotPlot-style visualization (uses AddModuleScore outputs) ----------
pathways_to_plot <- c(
  "apoptosis_pathway"
)
pathways_to_plot <- intersect(pathways_to_plot, colnames(obj@meta.data))
stopifnot(length(pathways_to_plot) == 1)

# summarize: mean score (color) & % cells with score > 0 (size)
dot_df <- obj@meta.data %>%
  dplyr::select(Samples, all_of(pathways_to_plot)) %>%
  tidyr::pivot_longer(all_of(pathways_to_plot), names_to = "Pathway", values_to = "Score") %>%
  dplyr::group_by(Samples, Pathway) %>%
  dplyr::summarise(mean_score = mean(Score, na.rm = TRUE),
                   pct_pos    = mean(Score > 0, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::group_by(Pathway) %>%
  dplyr::mutate(z_mean = scale(mean_score)[,1]) %>%
  dplyr::ungroup()

dot_df$Samples <- factor(dot_df$Samples, levels = sample_order)
dot_df$Pathway <- factor(dot_df$Pathway, levels = pathways_to_plot)

ggplot(dot_df, aes(x = Samples, y = Pathway)) +
  geom_point(aes(size = pct_pos, color = z_mean)) +
  scale_size(range = c(2.5, 9), name = "Pct > 0") +
  scale_color_viridis_c(name = "Mean score (z)", option = "D") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Module scores (AddModuleScore) — DotPlot-style by sample", x = NULL, y = NULL)+coord_flip()


