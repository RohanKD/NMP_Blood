library(corrplot)

# obtain expression values of the genes I want to test rom the Seurat object
expression_matrix <- FetchData(PM, vars = c("tal1", "etv2", "tcf15", "meox1", "id3", "id1", "sox2", "msgn1"))
expression_matrix <- round(expression_matrix, 1)

cor_expression_matrix <- cor(expression_matrix)
corrplot(cor_expression_matrix, is.corr = T, method = "square")
# 
# 
# 
# # Replace 'cell_name' with the name of the cell you want to examine
# cell_name <- "AAACCCAAGAAACTAC-1"
# 
# # Replace 'gene_name' with the name of the gene you are interested in
# gene_name <- "sox2"
# 
# # Access the gene expression matrix from the Seurat object
# gene_expression_matrix <- GetAssayData(pbmc)
# 
# # Find the expression value of the specific gene in the specified cell
# expression_value <- gene_expression_matrix[gene_name, cell_name]
# 
# # Print the result
# print(expression_value)