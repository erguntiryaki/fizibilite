# Data used in following analysis can be retrieved from GEO with accession number: GSE100501
# Associated publication:
# Peterson, V.M., et al., Multiplexed quantification of proteins and transcripts in single cells. Nat Biotechnol, 2017. 35(10): p. 936-939.

# This analysis involves simultaneous mRNA and protein measurements from PBMC samples of 
# donor 2 and donor 3 in aferomentioned public dataset.

################### SET UP SEURAT OBJECT ################################
library(Seurat)

# Set up Seurat object and add ADT assay for donor 2
rna2 <- as.sparse(read.csv(file = "./data/reap/2/mrna_matrix.txt", sep = "\t"))
rna2 <- CreateSeuratObject(rna2, project = "donor2")

adt2 <- as.sparse(read.csv(file = "./data/reap/2/protein_matrix.txt", sep = "\t"))
adt2 <- CreateAssayObject(adt2)

rna2[["ADT"]] <- adt2

# Set up Seurat object and add ADT assay for donor 3
rna3 <- as.sparse(read.csv(file = "./data/reap/3/mrna_matrix.txt", sep = "\t"))
rna3 <- CreateSeuratObject(rna3, project = "donor3")

adt3 <- as.sparse(read.csv(file = "./data/reap/3/protein_matrix.txt", sep = "\t"))
adt3 <- CreateAssayObject(adt3)

rna3[["ADT"]] <- adt3
gc()

# Merge two Seurat objects
pbmc <- merge(rna2, y = rna3, add.cell.ids = c("donor2", "donor3"), project = "pbmc")

rm(rna2, rna3, adt2, adt3)
gc()

pbmc
## An object of class Seurat 
## 32786 features across 7488 samples within 2 assays 
## Active assay: RNA (32738 features, 0 variable features)
## 1 other assay present: ADT

Assays(pbmc)
## "RNA" "ADT"
DefaultAssay(pbmc)
## "RNA"

#%%%%%%%% 1. QC
#install.packages("remotes")
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
library(DoubletFinder)

# Check metadata
head(pbmc@meta.data, 5)

pbmc <- PercentageFeatureSet(pbmc, pattern = "MT-", col.name = "percent.mt", assay = "RNA")
head(pbmc@meta.data,5) # check metadata after new QC metric is added

pbmc <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]", col.name = "percent.ribo")

pbmc <- PercentageFeatureSet(pbmc, pattern = "^HB[^(P)]", col.name = "percent.hb")

pbmc <- PercentageFeatureSet(pbmc, pattern = "PECAM1|PF4", col.name = "percent.plat")

colnames(pbmc@meta.data)

## Visualize QC metrics
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", 
           "nCount_ADT", "nFeature_ADT")
VlnPlot(pbmc, features = feats, pt.size = 0.1, ncol = 3) + 
  NoLegend()

FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

# Filter low quality cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Standard analysis workflow
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- RunUMAP(pbmc, dims = 1:30)

#%%%%% Doublet detection via DoubletFinder
# define the expected number of doublet cellscells.
nExp <- round(ncol(pbmc) * 0.03)  # expect 3% doublets
pbmc <- doubletFinder_v3(pbmc, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

head(pbmc@meta.data,2)

# Plot predicted doublets on UMAP
DimPlot(pbmc, group.by = "DF.classifications_0.25_0.09_213")

# Filter the Seurat object such that only singlets remain
pbmc <- subset(pbmc, subset = DF.classifications_0.25_0.09_213 == "Singlet")
pbmc
## 32786 features across 6873 samples within 2 assays

###################### Process Protein Assay #############################

DefaultAssay(pbmc)
## RNA

# Normalize ADT Data
pbmc <- NormalizeData(pbmc, normalization.method = "CLR", margin = 2, assay = "ADT")

# Check Assay Keys
pbmc[["RNA"]]@key
## "rna_"
pbmc[["ADT"]]@key
## "adt_"

DimPlot(pbmc, group.by = "orig.ident")
pbmc <- RunTSNE(pbmc, dims = 1:30)
DimPlot(pbmc, reduction = "tsne")
pbmc <- pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc)
DimPlot(pbmc, reduction = "umap")

###################### scPred ##################
## install.packages("cli")
## install.packages("devtools")
## if (!requireNamespace("BiocManager", quietly = TRUE))
##   install.packages("BiocManager")

## BiocManager::install("SingleCellExperiment")
## devtools::install_github("immunogenomics/harmony")
## devtools::install_github("powellgenomicslab/scPred")


library(dplyr)
library(ggplot2)
library(scPred)

reference <- scPred::pbmc_1

reference
## 32838 features across 3500 samples within 1 assay

# Re run analysis pipeline for reference
reference <- reference %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% 
  RunPCA(verbose = F) %>% RunUMAP(dims = 1:30)

DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE) + NoAxes()
# Set the identity as louvain with resolution 0.3

pbmc <- SetIdent(pbmc, value = "RNA_snn_res.0.8")

pbmc <- pbmc%>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F) %>% 
  RunUMAP(dims = 1:30)

transfer.anchors <- FindTransferAnchors(reference = reference, query = pbmc, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = reference$cell_type, 
                            dims = 1:30)
pbmc <- AddMetaData(object = pbmc, metadata = predictions)

DimPlot(pbmc, group.by = "predicted.id", label = T, repel = T) + NoAxes()

ggplot(pbmc@meta.data, aes(x = RNA_snn_res.0.8, fill = predicted.id)) + geom_bar() + 
  theme_classic()

############## Visualization of Multimodal Analysis
rownames(pbmc[["ADT"]]) # check available antibodies and respective names

#### Visualize the expression of some canonical markers

# CD56 for NK cells
p1 <- FeaturePlot(pbmc, 'rna_NCAM1')
p2 <- FeaturePlot(pbmc, 'adt_CD56-TACATAAG')
p1 | p2

# NK
FeaturePlot(pbmc, features= c('NKG7', 'GNLY'))

# B cells
FeaturePlot(pbmc, c('CD79A', 'MS4A1'))
p1 <- FeaturePlot(pbmc, 'rna_MS4A1')
p2 <- FeaturePlot(pbmc, 'adt_CD20-ACGCGGAA')
p1 | p2

## T Cells
FeaturePlot(pbmc, c('CD3D', 'CD3G', 'CD3E'))
p1 <- FeaturePlot(pbmc, 'rna_CD3D')
p2 <- FeaturePlot(pbmc, 'adt_CD3-AGGATCGA')
p1 | p2

# CD4 T lymphocytes
FeaturePlot(pbmc, c('IL7R', 'CCR7'))
p1 <- FeaturePlot(pbmc, 'rna_CD4')
p2 <- FeaturePlot(pbmc, 'adt_CD4-CACGATTC')
p1 | p2

# CD8 T lymphocytes
p1 <- FeaturePlot(pbmc, 'rna_CD8A')
p2 <- FeaturePlot(pbmc, 'adt_CD8a-ACCCGCAC')
p1 | p2
