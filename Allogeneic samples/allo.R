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
  Control = "D:/new data/xilis/h1 colon and t cell single cell RNA seq/colon/raw_matrix_h1_colon_control",
  Utr_T = "D:/new data/xilis/h1 colon and t cell single cell RNA seq/colon/colon_utr_raw_matrix",
  Her2_T = "D:/new data/xilis/h1 colon and t cell single cell RNA seq/colon/colon_her2_raw_matrix"
)
seurat_list <- list()


# Load, QC, Normalize
for (sample_name in names(sample_paths)) {
  cat("Processing:", sample_name, "\n")
  
  # Read 10X
  counts <- Read10X(data.dir = sample_paths[[sample_name]])
  obj <- CreateSeuratObject(counts = counts, project = sample_name)
  obj$orig.ident <- sample_name
  
  # QC
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > 250 & nFeature_RNA < 7500 & percent.mt < 4)
  
  # Normalize, variable genes, PCA
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20)
  
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:20)
  obj <- FindClusters(obj, resolution = 1)
  
  seurat_list[[sample_name]] <- obj
  
}

# Merge all samples into one Seurat object
Colon_T<- merge(x =seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "Merged")


# Normalize, variable genes, scale, PCA again on merged object
Colon_T <- NormalizeData(Colon_T)
Colon_T <- FindVariableFeatures(Colon_T, selection.method = "vst", nfeatures = 3000)
Colon_T <- ScaleData(Colon_T, verbose = FALSE)

Colon_T <- RunPCA(Colon_T, npcs = 50, verbose = FALSE)
ElbowPlot(Colon_T, ndims = ncol(Embeddings(Colon_T, "pca")))
PCHeatmap(Colon_T, dims = 1:10, cells = 500, balanced = TRUE, ncol = 4)
Colon_T <- FindNeighbors(Colon_T, reduction = "pca", dims = 1:30)
Colon_T <- FindClusters(Colon_T, resolution = 1)
Colon_T <- RunUMAP(Colon_T, reduction = "pca", dims = 1:30)

harmony <- IntegrateLayers(
  object = Colon_T, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

harmony <- FindNeighbors(harmony, reduction = "harmony", dims = 1:30)
harmony <- FindClusters(harmony, resolution = 1.4, cluster.name = "harmony_clusters")
harmony <- RunUMAP(harmony, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(
  harmony,
  reduction = "umap",
  combine = F, label.size = 2,label = TRUE,split.by = "orig.ident")

DimPlot(
  harmony,
  reduction = "umap",
  combine = F, label.size = 2,label = TRUE)



harmony<- JoinLayers(harmony)
harmony_markers <- FindAllMarkers(harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
harmony_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10_harmony_markers <- harmony_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(harmony, features = top10_harmony_markers$gene) + NoLegend()
write.csv(harmony_markers, "D:/new data/xilis/h1 colon and t cell single cell RNA seq/colon/harmony_markers_new.csv",row.names = T)

table(harmony$orig.ident)

# cell anotation 
harmony_ident <- setNames(c("Enterocyte","Enterocyte","Enterocyte","Goblet cell","Enterocyte",
                            "Unknown","T cell","T cell","T cell","Fibroblast","Fibroblast",
                            "T cell","Enterocyte","Enterocyte","Enteroendocrine cell","Enteroendocrine cell",
                            "Proliferating cell","Neuronal cell","Neuronal cell"),
                          levels(harmony))

harmony<- RenameIdents(harmony, harmony_ident)
harmony$clusters <- Idents(harmony)

DimPlot(
  harmony,
  reduction = "umap",
  group.by = c("clusters"),
  combine = F, label.size = 2,split.by="orig.ident"
)

harmony_control <- subset(harmony, subset = orig.ident == "Control")
harmony_utr <- subset(harmony, subset = orig.ident == "Utr_T")
harmony_her2 <- subset(harmony, subset = orig.ident == "Her2_T")
table(harmony_control[["clusters"]])


harmony <- FindClusters(harmony, resolution = 1.4, cluster.name = "harmony_clusters")
harmony <- RunUMAP(harmony, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
harmony_ident1 <- setNames(c("Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop",
                             "Colon_iDrop","T_cell","T_cell","T_cell","Colon_iDrop","Colon_iDrop",
                             "T_cell","Colon_iDrop","Colon_iDrop","Colon_iDrop","Colon_iDrop",
                             "Colon_iDrop","Colon_iDrop","Colon_iDrop"),
                          levels(harmony))
harmony<- RenameIdents(harmony, harmony_ident1)
harmony$Immu_iDrop <- Idents(harmony)

DimPlot(harmony, reduction = "umap", label = T)+NoLegend()

DimPlot(harmony, reduction = "umap", label = TRUE,split.by = "orig.ident",repel = T)

# Check
table(harmony$orig.ident)

FeaturePlot(harmony, features = c("CD3D", "EPCAM", "DCN","PAX6"), reduction = "umap") & 
  xlab("UMAP_1") &
  ylab("UMAP_2") & 
  theme(axis.title = element_text(size = 14))

VlnPlot(harmony,c("CD3D","CD3E","EPCAM","CDH1","DCN","COL1A1","PAX6","NEFL"),ncol = 4)

#samples

Idents(harmony) <- "orig.ident"
harmony$Samples <- Idents(harmony)

sample_colors <- c(
  "Control" = "#0a9396",
  "Utr_T" = "#ee9b00",
  "Her2_T" = "#ae2012"
)


DimPlot(
  harmony,
  group.by = "Samples",
  label = F,
  cols = sample_colors,
  pt.size = 0.5
)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 14, face = "bold"))

#two clusters

Idents(harmony) <- "Immu_iDrop"

DimPlot(
  harmony,
  group.by = "Immu_iDrop",
  label = T,
  pt.size = 0.5
)+NoLegend()+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 14, face = "bold"))


#T cells
harmony_T <- subset(harmony, subset = Immu_iDrop == "T_cell")

Idents(harmony_T) <- "Samples"

de_Her2_vs_Untr_T <- FindMarkers(
  harmony_T,
  ident.1 = "Her2_T",
  ident.2 = "Utr_T",
  logfc.threshold = 0.25,
  test.use = "wilcox"  # default
)

sig_de_Her2_vs_Untr_T <- de_Her2_vs_Untr_T %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

head(sig_de_Her2_vs_Untr_T)

write.csv(sig_de_Her2_vs_Untr_T, "D:/new data/xilis/h1 colon and t cell single cell RNA seq/colon/sig_de_Her2_vs_Untr_T.csv",row.names = T)


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

table(harmony_T$orig.ident)

harmony_T <- subset(harmony_T, subset = orig.ident %in% c("Utr_T", "Her2_T"))

DotPlot(harmony_T, c("TNF","TNFSF10","LTA","GZMB","CSF1","CSF2","IL13","IL5",
   "IFNG","EIF2AK2","ICAM1","ISG15","OAS1","CXCL10","CXCL9",
   "CCR4","CCR7","CXCL11","CXCL5","GNG8","XCL1"))+
   theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y  = element_text(face = "italic"))+coord_flip()

DotPlot(harmony_T, c("FAS","TNFRSF10A","TNFRSF10B","CASP3","CASP1","IRF1",
                     "IRF3"),group.by = "Samples")+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y  = element_text(face = "italic")) +coord_flip()

#colon iDrop
harmony_Colon <- subset(
  harmony,
  idents = c("T_cell"),
  invert = TRUE
)

Idents(harmony_Colon) <- "clusters"

DimPlot(harmony_Colon, pt.size = 0.5)+
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
harmony_Colon$samples <- factor(harmony_Colon$orig.ident, levels = c("Control", "Utr_T", "Her2_T"))

# Assign to cluster identities
names(colors16) <- levels(harmony_Colon)


DimPlot(harmony_Colon,cols = colors16,pt.size = 0.5)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 14))

DimPlot(harmony_Colon, pt.size = 0.5,split.by = "Samples",cols = colors16)+
  xlab("UMAP_1") + 
  ylab("UMAP_2") +
  theme(axis.title = element_text(size = 14))


# subset
harmony_Colon_control <- subset(harmony_Colon, subset = orig.ident == "Control")

table(harmony_Colon_control[["clusters"]])

# Extract ERBB2 expression and metadata
plot_data <- FetchData(harmony_Colon_control, vars = c("ERBB2", "clusters", "orig.ident"))

# Filter out ERBB2 = 0
plot_data <- plot_data %>% filter(ERBB2 > 0)
table(plot_data[["clusters"]])

# Ensure 'celltype' is a factor with correct ordering
plot_data$celltype <- factor(plot_data$clusters)

# Create named color vector for all cell types
n_ct <- length(unique(plot_data$celltype))

# Either use predefined palette if enough colors
if (n_ct <= 16) {
  my_colors <- c(
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
  )[1:n_ct]
} else {
  my_colors <- scales::hue_pal()(n_ct)  # auto-generate if >16
}

celltype_colors <- setNames(my_colors, levels(plot_data$celltype))

# Create the jitter plot
ggplot(plot_data, aes(x = celltype, y = ERBB2, color = celltype)) +
  geom_boxplot()+
  geom_jitter(width = 0.25, size = 0.6, alpha = 0.6) +
  scale_color_manual(values = celltype_colors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = c("bold.italic"))
  ) +
  labs(
    title = "ERBB2",
    x = "Cell Type",
    y = "Normalized Expression"
  )

# her2 vs control

Idents(harmony_Colon) <- harmony_Colon$orig.ident
table(Idents(harmony_Colon))  # Check counts


de_Her2_vs_Control_harmony_Colon <- FindMarkers(
  harmony_Colon,
  ident.1 = "Her2_T",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  test.use = "wilcox"  # default
)

sig_de_Her2_vs_Control_harmony_Colon <- de_Her2_vs_Control_harmony_Colon %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

head(sig_de_Her2_vs_Control_harmony_Colon)

write.csv(sig_de_Her2_vs_Control_harmony_Colon, "D:/new data/xilis/h1 colon and t cell single cell RNA seq/colon/sig_de_Her2_vs_Control_harmony_Colon.csv",row.names = T)


de_Her2_vs_Control_harmony_Colon$gene <- rownames(de_Her2_vs_Control_harmony_Colon)

de_Her2_vs_Control_harmony_Colon$sig <- with(de_Her2_vs_Control_harmony_Colon,
                                             ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, "Significant", "Not"))

# Example: adjust this to your DE result
res_df <- de_Her2_vs_Control_harmony_Colon
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
  slice(1:10)

# Plot
ggplot(res_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = sig), alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = gene),size=4,max.overlaps = Inf,fontface = "italic") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Colon_iDrop: Her2_T vs Control",
    x = "logFC",
    y = "-log10(p)"
  ) +
  theme(legend.position = "top")

#utr vs control
Idents(harmony_Colon) <- harmony_Colon$orig.ident
table(Idents(harmony_Colon))  # Check counts


de_Utr_vs_Control_harmony_Colon <- FindMarkers(
  harmony_Colon,
  ident.1 = "Utr_T",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  test.use = "wilcox"  # default
)

sig_de_Utr_vs_Control_harmony_Colon <- de_Utr_vs_Control_harmony_Colon %>%
  dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

head(sig_de_Utr_vs_Control_harmony_Colon)

write.csv(sig_de_Utr_vs_Control_harmony_Colon, "D:/new data/xilis/h1 colon and t cell single cell RNA seq/colon/sig_de_Utr_vs_Control_harmony_Colon.csv",row.names = T)


de_Utr_vs_Control_harmony_Colon$gene <- rownames(de_Utr_vs_Control_harmony_Colon)

de_Utr_vs_Control_harmony_Colon$sig <- with(de_Utr_vs_Control_harmony_Colon,
                                            ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, "Significant", "Not"))

# Example: adjust this to your DE result
res_df <- de_Utr_vs_Control_harmony_Colon
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
  geom_text_repel(data = top_genes, aes(label = gene),size=4,max.overlaps = Inf,fontface = "italic") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Colon_iDrop: Utr_T vs Control",
    x = "logFC",
    y = "-log10(p)"
  ) +
  theme(legend.position = "top")



DotPlot(harmony_Colon, c("BCL2","RIPK1","BAK1","BBC3","NFKBIA",
                         "PMAIP1","FAS","FADD","TNFRSF10B",
                         "TNFRSF1B","CASP8",
                         "CASP7","CASP4","CASP1","IRF7",
                         "IRF3","IRF2","IRF1","BIRC3"),group.by = "Samples")+
  theme(axis.text.x=element_text(angle=45,hjust = 1),
        axis.text.y = element_text(face = "italic"))+coord_flip()









