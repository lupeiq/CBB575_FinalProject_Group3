# CBB575_FinalProject_Group3
# Single-Cell RNA Sequencing Analysis of Influenza-Infected Lung Tissue

## Overview
This project investigates the cellular and molecular dynamics of lung tissue in response to influenza infection using single-cell RNA sequencing (scRNA-seq). By analyzing control and influenza-infected mouse lungs, we identify key immune pathways, transcriptional changes, and cell-type interactions to gain insights into the immune response and disease progression.

## Features
- **Data Preprocessing**: Quality control, normalization, and dimensionality reduction using the Seurat R package.
- **Cell Clustering and Annotation**: Identification of cell populations with UMAP visualizations and cell-type annotation through transfer learning.
- **Differential Gene Expression Analysis**: Detection of differentially expressed genes (DEGs) using the MAST method.
- **Pathway Enrichment**: Exploration of immune and non-immune pathways via pathway analysis.
- **Immune Cell Focus**: Detailed analysis of NK cells and CD4 T cells, including cluster-specific DEGs and pathway interactions.

## Workflow
1. **Raw Data Preprocessing**:
   - Perform scRNA-seq data quality control using Cell Ranger and Seurat.
   - Filter cells based on mitochondrial gene percentage and feature count thresholds.
   
2. **Clustering and Dimensional Reduction**:
   - Use PCA for linear dimensionality reduction.
   - Perform clustering with the Louvain algorithm and visualize with UMAP.

3. **Annotation**:
   - Annotate cell types by transferring metadata from a reference dataset.
   - Visualize annotated clusters using customized UMAP plots.

4. **Differential Expression Analysis**:
   - Identify DEGs between control and infected samples using the MAST framework.
   - Filter DEGs with |log2FC| > 1 and adjusted p-value < 0.05.
   - Create volcano plots for visualization.

5. **Immune Cell Subpopulation Analysis**:
   - Focus on NK cells and CD4 T cells for detailed cluster-specific analysis.
   - Investigate pathways enriched in different cell clusters.

6. **Pathway Enrichment Analysis**:
   - Use MetaCore for pathway enrichment analysis of identified DEGs.
   - Highlight immune-related and non-immune pathways (e.g., neuronal signaling).

## Key Findings
- Increased immune cell populations, including CD4 T, CD8 T, and NK cells, in infected lung tissue.
- Identification of 2,003 DEGs, including immune-related pathways and interactions.
- Novel findings include NK cell interactions with epithelial cells and CD4 T cell involvement in neuronal pathways.

## Requirements
### Software
- [Cell Ranger 7.1.0](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome)
- R (≥ 4.0.0) with the following libraries:
  - `Seurat`
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `MAST`

### Data
- Available in [Google Drive](https://drive.google.com/drive/folders/1lOUGvbXTEetbv05i1dEYK0U1QdEVnTg3?usp=drive_link)

## Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/lupeiq/CBB575_FinalProject_Group3.git
   cd CBB575_FinalProject_Group3

