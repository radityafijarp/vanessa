#install.packages("pheatmap")
library(pheatmap)
pheatmap(as.matrix(iris.data), 
         cluster_rows = iris.tree,
         cluster_cols = TRUE, 
         scale = "column", 
         labels_row = as.character(iris.species),
         fontsize_row = 3)

# Create an annotation data.frame with the same row names as 
# the input data (iris.data)
rownames(iris.data) <- 1:150 # OBS: iris did not have row names
annot.row <- data.frame(row.names = rownames(iris.data), Species = iris.species)

pheatmap(as.matrix(iris.data), 
         cluster_rows = iris.tree,
         cluster_cols = TRUE, 
         scale = "column", 
         show_rownames = FALSE, 
         annotation_row = annot.row)

# Transpose the data to cluster samples
aspwood.dist <- dist(t(aspwood), method = "euclidean")
aspwood.tree <- hclust(aspwood.dist, method = "complete")

# Create an annotation data frame for the stages
sample_names <- colnames(aspwood)
annot.col <- data.frame(row.names = sample_names, Stage = aspwood.stages)

# Plot the heatmap
pheatmap(as.matrix(aspwood), 
         cluster_rows = TRUE, # Default clustering for rows (genes)
         cluster_cols = aspwood.tree, # Precomputed clustering for columns (samples)
         scale = "row", # Scale data across rows (genes)
         annotation_col = annot.col, # Color samples by stage
         show_rownames = FALSE, # Hide row names for readability
         show_colnames = TRUE) # Show column names (sample names)


# Plot the heatmap
pheatmap(as.matrix(aspwood), 
         cluster_rows = TRUE, # Default clustering for rows (genes)
         cluster_cols = aspwood.tree, # Precomputed clustering for columns (samples)
         annotation_col = annot.col, # Color samples by stage
         show_rownames = FALSE, # Hide row names for readability
         show_colnames = TRUE) # Show column names (sample names)

# Effect of Removing scale = "row":

# If you omit scale = "row", the heatmap will represent raw expression values.
# High-expressing genes will dominate the visualization, making it hard to
# compare gene expression patterns across samples. Including scale = "row"
# allows easier identification of up/down-regulated genes across samples