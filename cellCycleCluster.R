library(AnnotationHub)
library(dplyr)
cell_cycle_genes <- read.csv("C:\\Users\\rohan\\OneDrive\\Documents\\gene_analysis\\Danio_rerio.csv")
ah <- AnnotationHub()
ahDb <- query(ah, 
              pattern = c("Danio rerio", "EnsDb"), 
              ignore.case = TRUE)
# if (!require ("BiocManager", quietly = TRUE)) install.packages ("BiocManager")
# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring
pbmc <- CellCycleScoring(pbmc,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)

# Perform PCA and color by cell cycle phase
pbmc <- RunPCA(pbmc
               )

# Visualize the PCA, grouping by cell cycle phase
DimPlot(pbmc,
        reduction = "pca",
        group.by= "Phase")

# alternater G2/S regression

marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score
marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))