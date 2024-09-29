
# Load the libraries
library(tidyverse)
library(usedist)
library(dendextend)
library(pheatmap)

# Save our own version of the iris data
iris.data <- iris

# Save the species data
iris.species <- iris.data$Species

# Remove the species information so that it doesn't influence the clustering
iris.data$Species <- NULL

library(usedist)

# Calculate the distance matrix 
# and show the first 5 rows/columns
iris.data %>% 
  dist() %>% 
  dist_subset(1:5)

# What distance function does dist use as default?
#  The dist function uses the Euclidean distance as its default method.

# Does the dist function calculate distances between rows or columns?
# The dist function calculates distances between rows.

# What hierarchical clustering method does hclust use as default?
# The hclust function uses the "complete" linkage method by default.
# Perform the hierarchical clustering
iris.tree <- iris.data %>% 
  dist() %>% 
  hclust()
# Plot the similarity tree (dendrogram)
iris.tree %>% 
  as.dendrogram() %>% 
  set("labels", iris.species[iris.tree$order]) %>% # change labels
  set("labels_cex", 0.25) %>% # label size
  set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 0.5) %>%  # node point size
  set("leaves_col", iris.species[iris.tree$order]) %>%  # node point color
  set("branches_k_color", k=3) %>% # branch color
  plot()

# Difference Between Coloring of Branches and Leaves
# The branches are colored based on the three main clusters identified by the hierarchical clustering algorithm, while the leaves are colored according to the actual species labels of the iris flowers. This allows us to compare how well the clustering algorithm grouped the data in comparison to the true species classifications.


# Annotate the stages of wood formation
aspwood.stages <- factor(c(rep("PhloemZone", 5), 
                           rep("ExpansionZone", 5), 
                           rep("SecondaryCellWallZone", 9), 
                           rep("LignificationZone", 6)))

# Calculate the distance matrix using the transposed data
aspwood.dist <- aspwood %>% 
  t() %>% 
  dist()

# Perform hierarchical clustering
aspwood.tree <- hclust(aspwood.dist)

# Reorder the tree to match the sampling order (1 to 25)
aspwood.tree <- aspwood.tree %>% 
  as.dendrogram() %>% 
  reorder(1:25, agglo.FUN = mean) %>% 
  as.hclust()

# Plot the hierarchical clustering dendrogram
aspwood.tree %>% 
  as.dendrogram() %>% 
  set("labels", aspwood.stages[aspwood.tree$order]) %>% # Change labels to stages
  set("labels_cex", 0.5) %>% # Adjust label size
  set("leaves_pch", 19) %>%  # Node point type
  set("leaves_cex", 0.7) %>%  # Node point size
  set("leaves_col", aspwood.stages[aspwood.tree$order]) %>% # Node point color
  set("branches_k_color", k = 4) %>% # Branch color (you can change k to reflect the number of clusters)
  plot()

# Biological Interpretation of the Tree
# Clusters: Observe how samples from different stages group together. Samples from the same stage (e.g., PhloemZone, ExpansionZone) should ideally cluster together, indicating that their gene expression profiles are similar.

# Branching: The tree's structure shows how closely related different stages are in terms of their gene expression. Shorter branch lengths indicate higher similarity.

# Comparison with the AspWood Paper: If the clusters closely match the paper's findings, it validates that the normalized expression data accurately capture the different stages of wood formation.