library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:\\Users\\rohan\\OneDrive\\Documents\\gene_analysis\\full_tailbud")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

#Identifying "marker genes" top 2k
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)


# scale data: The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
# 
#Principal component analysis!!!
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# visualize PCA
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# DimPlot(pbmc, reduction = "pca")
#Heatmap yellow/purple
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# elbow plot is so much better
# ElbowPlot(pbmc)
# will use 18 PCs
# 
# clustering algorithm
pbmc <- FindNeighbors(pbmc, dims = 1:17)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
# head(Idents(pbmc), 5)
# 
#Clustering plot UMAP
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:17)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
# 
# #Identify each cluster using marker identification and biological knowledge
# # find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)
# 
# # find all biomarkers
# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
# 
# # Seurat has several tests for differential expression which can be set with the test.use parameter 
# # (see our DE vignette for details). 
# # For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
# 
# cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# 
# # plot violin plot to show expression probability distributions across clusters
# # alter feature (gene) choices - maybe try tbxta + sox2
# VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# 
# # you can plot raw counts as well
# VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# 
# #Shows gene expressionacross cell population for different genes (same as David's plots)
# # switch out gene choices
FeaturePlot(pbmc, features = c("sox2"))
# 
# # Plot markers and generate marker gene heatmap for top 20 genes for each cluster
pbmc.markers %>%
  # group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
# 
# 
# # find markers for each cluster and try to manually label, use david marker gene info and label 15 clusters
# # Implement looking for markers that are coexpressing paraxial mesoderm and endothelial + blood precursor
PM_filtered <- subset(PMf, subset = meox1 > 1 & myf5>1 & tal1 >1)
# # find what other genes are expressed
# 
# # run this code for every cluster
# cluster2.markers <- FindMarkers(pbmc, ident.1 = change_this, min.pct = 0.25, max.cells.per.ident = 1000)
# head(cluster2.markers, n = 10)
new.cluster.ids <- c("0", "1", "Blood related", "Mesodermal progenitors", "Neuron related", "NT1",
                     "Anterior PM", "7", "8", "Notochord", "10", "NMPs", "blood related", "epidermis", "IM", "NT2", "Periderm1","Periderm2", "Myoblasts")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# Step 1: Create a list of gene names you want to filter by
gene_list <- c("etv2", "meox1", "myf5")

# Step 2: Calculate the expression values for the genes in your Seurat object
# Replace 'your_data' with the actual gene expression data in your Seurat object
gene_expression <- FetchData(PMf, vars = gene_list)

# Step 3: Use the Subset function to filter the cells based on gene expression values
blood_precursor <- subset(PMf, subset = gene_expression$gene1 > 1 & 
                               gene_expression$gene2 > 1 & 
                               gene_expression$gene3 > 1)