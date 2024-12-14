########## raw data creation (We will use the data after this step.) ########## 

#!/bin/bash

#SBATCH --time=24:00:00 -p day --ntasks=1 --cpus-per-task=8 --mem=128G --job-name=Sample1.Cellranger -o /vast/palmer/scratch/xiting_yan/pc775/Cellranger/script/Sample1_cellranger.sh.o%J -e /vast/palmer/scratch/xiting_yan/pc775/Cellranger/script/Sample1_cellranger.sh.e%J

#ml CellRanger/7.1.0
#cd /vast/palmer/scratch/xiting_yan/pc775/Cellranger/sample_out

#cellranger count --id=Sample1 --sample=Sample1 --transcriptome=/gpfs/gibbs/pi/xiting_yan/pc775/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0 --fastqs=/vast/palmer/scratch/xiting_yan/pc775/multiOme/sample_out/Sample1/RNA --chemistry=ARC-v1 --localcores=20 --localmem=128



########## Seurat analysis ########## 
### load library
library(dplyr)
library(tidyr)
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
library(Seurat)
library(patchwork)
library(ggplot2)

set.seed(1234)
### Load scRNAseq data
RNA1.data <- Read10X(data.dir = "/home/cbb575_pc775/final_project/s1_filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
RNA1 <- CreateSeuratObject(counts = RNA1.data, project = "RNA1", min.cells = 3, min.features = 200)

### QC and selecting cells for further analysis (The [[ operator can add columns to object metadata.)
RNA1[["percent.mt"]] <- PercentageFeatureSet(RNA1, pattern = "mt-")
VlnPlot(RNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RNA1 <- subset(RNA1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

### Normalizing data
RNA1 <- NormalizeData(RNA1, normalization.method = "LogNormalize", scale.factor = 10000)
RNA1 <- NormalizeData(RNA1)

### Scaling the data
all.genes <- rownames(RNA1)
RNA1 <- ScaleData(RNA1, features = all.genes)

### Perform linear dimensional reduction
RNA1 <- FindVariableFeatures(object = RNA1)
RNA1 <- RunPCA(RNA1, features = VariableFeatures(object = RNA1))

### Determine the ‘dimensionality’ of the dataset
ElbowPlot(RNA1)

### Cluster the cell
RNA1 <- FindNeighbors(RNA1, dims = 1:20)
RNA1 <- FindClusters(RNA1, resolution = 0.5)

### Run non-linear dimensional reduction (UMAP/tSNE)
RNA1 <- RunUMAP(RNA1, dims = 1:20)
saveRDS(RNA1, file = "/home/cbb575_pc775/final_project/s1_SeuratObject.rds")



RNA3.data <- Read10X(data.dir = "/home/cbb575_pc775/final_project/s3_filtered_feature_bc_matrix")
RNA3 <- CreateSeuratObject(counts = RNA3.data, project = "RNA3", min.cells = 3, min.features = 200)
RNA3[["percent.mt"]] <- PercentageFeatureSet(RNA3, pattern = "mt-")
VlnPlot(RNA3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
RNA3 <- subset(RNA3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
RNA3 <- NormalizeData(RNA3, normalization.method = "LogNormalize", scale.factor = 10000)
RNA3 <- NormalizeData(RNA3)
all.genes <- rownames(RNA3)
RNA3 <- ScaleData(RNA3, features = all.genes)
RNA3 <- FindVariableFeatures(object = RNA3)
RNA3 <- RunPCA(RNA3, features = VariableFeatures(object = RNA3))
ElbowPlot(RNA3)
RNA3 <- FindNeighbors(RNA3, dims = 1:20)
RNA3 <- FindClusters(RNA3, resolution = 0.5)
RNA3 <- RunUMAP(RNA3, dims = 1:20)
saveRDS(RNA3, file = "/home/cbb575_pc775/final_project/s3_SeuratObject.rds")



########## Annotation ##########
RNA.ref <- readRDS("/home/cbb575_pc775/final_project/mouselung.RDS")

### pre-process dataset (without integration)
RNA.ref <- NormalizeData(RNA.ref)
RNA.ref <- FindVariableFeatures(RNA.ref)
RNA.ref <- ScaleData(RNA.ref)
RNA.ref <- RunPCA(RNA.ref)
RNA.ref <- FindNeighbors(RNA.ref, dims = 1:20)
RNA.ref <- FindClusters(RNA.ref)
RNA.ref <- RunUMAP(RNA.ref, dims = 1:20)
DimPlot(RNA.ref, group.by = c("celltype_level3"))

### load query data
RNA.query <- readRDS("/home/cbb575_pc775/final_project/s1_SeuratObject.rds")
# RNA.query <- RNA1
RNA.query <- NormalizeData(RNA.query)
RNA.anchors <- FindTransferAnchors(reference = RNA.ref, query = RNA.query, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = RNA.anchors, refdata = RNA.ref$celltype_level3, dims = 1:20)
RNA.query <- AddMetaData(RNA.query, metadata = predictions)

### UMAP plots
my.distinct.colors = c(
  "#e6194b", "#3cb44b", "#4363d8", "#ffd8b1", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080", 
  "#9a6324", "#800000", "#aaffc3", "#808000", "#f58231", "#808080", "#e6beff", "#000075", "#ffe119", "#000000", 
  "#00ff00", "#ff4500", "#00ced1", "#556b2f", "#a0522d", "#8b0000", "#808000", "#483d8b", "#008000", "#008080",
  "#4682b4", "#000080", "#9acd32", "#daa520", "#7f007f", "#8fbc8f", "#b03060", "#d2b48c", "#696969", "#ff8c00", 
  "#00ff7f", "#dc143c", "#f4a460", "#0000ff", "#a020f0", "#adff2f", "#ff00ff", "#1e90ff", "#f0e68c", "#fa8072",
  "#ffff54", "#dda0dd", "#87ceeb", "#7b68ee", "#ee82ee", "#98fb98", "#7fffd4", "#ffb6c1", "#dcdcdc", "#000000",
  "#00ff00", "#ff4500", "#00ced1", "#556b2f", "#a0522d", "#8b0000", "#808000", "#483d8b", "#008000", "#008080",
  "#4682b4", "#000080", "#9acd32", "#daa520", "#7f007f", "#8fbc8f", "#b03060", "#d2b48c", "#696969", "#ff8c00", 
  "#00ff7f", "#dc143c", "#f4a460", "#0000ff", "#a020f0", "#adff2f", "#ff00ff", "#1e90ff", "#f0e68c", "#fa8072",
  "#ffff54", "#dda0dd", "#87ceeb", "#7b68ee", "#ee82ee", "#98fb98", "#7fffd4", "#ffb6c1", "#dcdcdc", "#000000")

DimPlotUMAP_name = "s1_DimPlotUMAP_annotation.png"
png(DimPlotUMAP_name, res=600, height=5000, width=9000)
DimPlot(RNA.query, group.by = "predicted.id", label = FALSE) +
    theme(
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.text = element_text(size = 10),  # Adjust legend text size
        legend.title = element_text(size = 12) # Adjust legend title size
    ) +
    guides(
        color = guide_legend(
            nrow = 3, 
            byrow = TRUE,
            override.aes = list(size = 3)  # Adjust circle size in legend
        )
    ) +
  scale_color_manual(values = my.distinct.colors)
dev.off()
# Save file
saveRDS(RNA.query, file = "/home/cbb575_pc775/final_project/s1_SeuratObject_annotation.rds")



RNA.query <- readRDS("/home/cbb575_pc775/final_project/s3_SeuratObject.rds")
# RNA.query <- RNA3
RNA.query <- NormalizeData(RNA.query)
RNA.anchors <- FindTransferAnchors(reference = RNA.ref, query = RNA.query, dims = 1:20, reference.reduction = "pca")
predictions <- TransferData(anchorset = RNA.anchors, refdata = RNA.ref$celltype_level3, dims = 1:20)
RNA.query <- AddMetaData(RNA.query, metadata = predictions)
DimPlotUMAP_name = "s3_DimPlotUMAP_annotation.png"
png(DimPlotUMAP_name, res=600, height=5000, width=9000)
DimPlot(RNA.query, group.by = "predicted.id", label = FALSE) +
    theme(
        legend.position = "bottom", 
        plot.title = element_blank(),
        legend.text = element_text(size = 10),  # Adjust legend text size
        legend.title = element_text(size = 12) # Adjust legend title size
    ) +
    guides(
        color = guide_legend(
            nrow = 3, 
            byrow = TRUE,
            override.aes = list(size = 3)  # Adjust circle size in legend
        )
    ) +
  scale_color_manual(values = my.distinct.colors)
dev.off()
saveRDS(RNA.query, file = "/home/cbb575_pc775/final_project/s3_SeuratObject_annotation.rds")



######### Merge sample1 and sample3 seurat objects #########
RNA1 <- readRDS("/home/cbb575_pc775/final_project/s1_SeuratObject_annotation.rds")
RNA3 <- readRDS("/home/cbb575_pc775/final_project/s3_SeuratObject_annotation.rds")
RNA.combined <- merge(RNA1, y = RNA3, add.cell.ids = c("s1", "s3"), project = "RNA", merge.data = TRUE)
RNA.combined <- NormalizeData(RNA.combined)
RNA.combined <- FindVariableFeatures(RNA.combined)
RNA.combined <- ScaleData(RNA.combined)
RNA.combined <- RunPCA(RNA.combined)
RNA.combined <- FindNeighbors(RNA.combined, dims = 1:20)
RNA.combined <- FindClusters(RNA.combined)
RNA.combined <- RunUMAP(RNA.combined, dims = 1:20)
saveRDS(RNA.combined, file = "/home/cbb575_pc775/final_project/combine_SeuratObject_annotation.rds")



########## Compare immune cells between two samples ##########
immune <- c("Basophil", "Neutrophil", "B", "CD4 T", "CD8 T", "Treg",
            "NK", "ILC", "Mast", "cDC1", "cDC2", "maDC")
lack.immune <- c("Basophil", "CD8 T", "Mast", "NK", "Treg")

RNA.immune_s1 <- RNA.combined@meta.data %>%
    dplyr::filter(predicted.id %in% immune) %>%
    dplyr::filter(orig.ident == "RNA1") %>%
    group_by(predicted.id) %>%
    summarise(cell_num = n()) %>%
    rename(cell_name = predicted.id) %>%
    complete(cell_name = lack.immune, fill = list(cell_num = 0)) %>%
    mutate(cell_proportion = (cell_num / sum(cell_num))*100) %>%
    mutate(cell_name = factor(cell_name, levels = c("Basophil", "Neutrophil", "B", "CD4 T", "CD8 T", 
                                                    "Treg", "NK", "ILC", "Mast", "cDC1", "cDC2", "maDC"))) %>%
    arrange(cell_name) %>%
    mutate(sample = "control")

RNA.immune_s3 <- RNA.combined@meta.data %>%
    dplyr::filter(predicted.id %in% immune) %>%
    dplyr::filter(orig.ident == "RNA3") %>%
    group_by(predicted.id) %>%
    summarise(cell_num = n()) %>%
    mutate(cell_proportion = (cell_num / sum(cell_num))*100) %>%
    rename(cell_name = predicted.id) %>%
    mutate(cell_name = factor(cell_name, levels = c("Basophil", "Neutrophil", "B", "CD4 T", "CD8 T", 
                                                    "Treg", "NK", "ILC", "Mast", "cDC1", "cDC2", "maDC"))) %>%
    arrange(cell_name) %>%
    mutate(sample = "flu")

RNA.immune <- rbind(RNA.immune_s1, RNA.immune_s3)

### Plot bar chart
my.distinct.colors <- c("#000000", "#808080", "#8b0000", "#ffd8b1",
                        "#46f0f0", "#9acd32", "#911eb4", "#008000",
                        "#0000ff", "#ffff54", "#ff8c00", "#e6194b")

barchart_name <- "Immune cell proportion across samples bar.png"
png(barchart_name, res = 300, height = 3000, width = 2000)
ggplot(RNA.immune, aes(x = sample, y = cell_proportion, fill = cell_name)) +
    geom_bar(stat = "identity") +
    labs(x = "Sample", y = "Cell Proportion (%)", fill = "Cell Type") +
    theme_minimal() +
    theme(
        panel.grid = element_blank(), # Remove the background grid
        axis.line = element_line(),  # Add x-axis and y-axis lines
        axis.text.x = element_text(hjust = 1)  # Optional rotation for x-axis labels
    ) +
    scale_fill_manual(values = my.distinct.colors) 
dev.off()

### Focus on CD4, CD8, Treg and NK (Since they are the biggest changes in these two samples)
focus.immune <- c("CD4 T", "CD8 T", "Treg", "NK")
RNA.focus.immune_s1 <- RNA.combined
cell.barcodes <- RNA.focus.immune_s1@meta.data %>% 
        dplyr::filter((predicted.id %in% focus.immune) & orig.ident == "RNA1") %>% 
        rownames %>%
        unlist
RNA.focus.immune_s1 <- subset(RNA.focus.immune_s1, cell = cell.barcodes)

RNA.focus.immune_s3 <- RNA.combined
cell.barcodes <- RNA.focus.immune_s3@meta.data %>% 
        dplyr::filter((predicted.id %in% focus.immune) & orig.ident == "RNA3") %>% 
        rownames %>%
        unlist
RNA.focus.immune_s3 <- subset(RNA.focus.immune_s3, cell = cell.barcodes)

### Plot UMAP
umap_name <- "s1_Focused Immune cell UMAP.png"
png(umap_name, res = 300, height = 2000, width = 3000)
DimPlot(RNA.focus.immune_s1, group.by = c("predicted.id"), label = FALSE) +
    ggtitle("") 
dev.off() # CD4 are gathered together

umap_name <- "s3_Focused Immune cell UMAP.png"
png(umap_name, res = 300, height = 2000, width = 3000)
DimPlot(RNA.focus.immune_s3, group.by = c("predicted.id"), label = FALSE) +
    ggtitle("") 
dev.off() # Two different groups in CD4 and NK



########## Find differentially expressed genes (DEGs) ##########
### Load library
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("MAST")
library(MAST)

Idents(RNA.combined) <- "orig.ident"
RNA.combined <- JoinLayers(RNA.combined)
# Positive markers in sample3
positive_markers <- FindMarkers(
    object = RNA.combined,
    ident.1 = "RNA3",
    ident.2 = "RNA1",
    only.pos = TRUE,
    min.pct = 0.1,
    test.use = "MAST")
positive_DEG <- positive_markers %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) 
# Output DEG
filename <- "positive_DEG_MAST.csv"
write.csv(positive_DEG, filename, row.names = TRUE)

### All differentially expressed markers
markers <- FindMarkers(
    object = RNA.combined,
    ident.1 = "RNA3",
    ident.2 = "RNA1",
    test.use = "MAST")
DEG <- markers %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1 & (pct.1 > 0.1 | pct.2 > 0.1)) 
# Save file
filename <- "gene_MAST.csv"
write.csv(markers, filename, row.names = TRUE)

filename <- "DEG_MAST.csv"
write.csv(DEG, filename, row.names = TRUE)

### Volcano plot
markers2 <- markers %>% filter(pct.1 > 0.1 | pct.2 > 0.1)
markers2$diffexpressed <- "NO"
markers2$diffexpressed[markers2$avg_log2FC > 1 & markers2$p_val_adj < 0.05] <- "UP"
markers2$diffexpressed[markers2$avg_log2FC < -1 & markers2$p_val_adj < 0.05] <- "DOWN"

p <- ggplot(data = markers2, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed)) +
     geom_vline(xintercept = c(-1, 1), col = "black", linetype = 'dashed') +
     geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') +
     geom_point(size = 1) +
     scale_color_manual(values = c("blue", "black", "red"), 
                        labels = c("Downregulated", "Not significant", "Upregulated")) +
     #coord_cartesian(ylim = c(0, 0.0001), xlim = c(-10, 10)) +
     labs(color = 'Regulation', #legend_title
          x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
     ggtitle("") +
     theme_minimal() +
     theme(panel.grid = element_blank(),
           axis.line = element_line())
# Export plots
volcano_name <- "VolcanoPlot.png"
png(volcano_name, res=300, height=1584, width=2073)
p
dev.off()



########## Focus on NK & CD4 T different clusters' DEGs) ##########
### NK
RNA.focus.immune_s3 <- RNA.combined
cell.barcodes <- RNA.focus.immune_s3@meta.data %>% 
        dplyr::filter((predicted.id %in% c("NK")) & orig.ident == "RNA3") %>% 
        rownames %>%
        unlist
RNA.focus.immune_s3 <- subset(RNA.focus.immune_s3, cell = cell.barcodes)

# Plot UMAP
umap_name <- "s3_NK UMAP.png"
png(umap_name, res = 300, height = 2000, width = 3000)
DimPlot(RNA.focus.immune_s3, group.by = c("predicted.id"), label = FALSE) +
    ggtitle("") 
dev.off() 

umap_name <- "s3_NK cluster UMAP.png"
png(umap_name, res = 300, height = 2000, width = 3000)
DimPlot(RNA.focus.immune_s3, group.by = c("seurat_clusters"), label = TRUE) +
    ggtitle("") 
dev.off() 

# Find DEGs
Idents(RNA.focus.immune_s3) <- "seurat_clusters"
cluster.markers <- FindMarkers(
    object = RNA.focus.immune_s3, 
    ident.1 = c(1,8,28,31), 
    ident.2 = c(6,9,10,11,14,15,16,20,23,25,27,29,30),
    min.pct = 0.1,
    test.use = "MAST")

cluster.DEG <- cluster.markers %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) 

# Save file
filename <- "NK cluster gene_MAST.csv"
write.csv(cluster.markers, filename, row.names = TRUE)

filename <- "NK cluster DEG_MAST.csv"
write.csv(cluster.DEG, filename, row.names = TRUE)

### CD4
RNA.focus.immune_s3 <- RNA.combined
cell.barcodes <- RNA.focus.immune_s3@meta.data %>% 
        dplyr::filter((predicted.id %in% c("CD4 T")) & orig.ident == "RNA3") %>% 
        rownames %>%
        unlist
RNA.focus.immune_s3 <- subset(RNA.focus.immune_s3, cell = cell.barcodes)

umap_name <- "s3_CD4 UMAP.png"
png(umap_name, res = 300, height = 2000, width = 3000)
DimPlot(RNA.focus.immune_s3, group.by = c("predicted.id"), label = FALSE) +
    ggtitle("") 
dev.off() 

umap_name <- "s3_CD4 cluster UMAP.png"
png(umap_name, res = 300, height = 2000, width = 3000)
DimPlot(RNA.focus.immune_s3, group.by = c("seurat_clusters"), label = TRUE) +
    ggtitle("") 
dev.off() 

Idents(RNA.focus.immune_s3) <- "seurat_clusters"
cluster.markers <- FindMarkers(
    object = RNA.focus.immune_s3, 
    ident.1 = c(10,15,30), 
    ident.2 = c(4,6,14,23,27),
    min.pct = 0.1,
    test.use = "MAST")
cluster.DEG <- cluster.markers %>% dplyr::filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1) 

filename <- "CD4 cluster gene_MAST.csv"
write.csv(cluster.markers, filename, row.names = TRUE)

filename <- "CD4 cluster DEG_MAST.csv"
write.csv(cluster.DEG, filename, row.names = TRUE)