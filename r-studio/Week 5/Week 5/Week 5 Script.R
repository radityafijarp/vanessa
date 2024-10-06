# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(matrixStats)

# Load the data
load(file = "TCGA.RData")

# Question: How many genes and samples are in the dataset?
# Answer: 
dim(counts)
# This will output the number of genes (rows) and samples (columns).

# Question: How many samples are there for each class?
# Answer: 
table(classes)
# This will show how many samples there are for each tissue class.

# Select classes 'Colon', 'Rectum', 'Stomach'
idx <- classes %in% c("Colon", "Rectum", "Stomach")
classes <- factor(classes[idx])
counts <- counts[, idx]

# Filter genes with at least 1024 mapped reads in at least 10 samples
counts <- counts[rowSums(counts > 2^10) > 10,]
dim(counts)
# Output will show the new dimensions of the filtered dataset.

# Normalize using DESeq2 to compute VST expression values
dds <- DESeqDataSetFromMatrix(countData = counts, colData = data.frame(classes), design = ~ classes)
vst <- vst(dds)
expr <- assay(vst)

# Calculate standard deviations with useNames = FALSE (no need to retain row names for now)
standard.deviations <- rowSds(expr, useNames = FALSE)


# Visualizing the distribution of standard deviations
tibble(StandardDeviation = standard.deviations) %>%
  ggplot(aes(x = StandardDeviation)) +
  geom_density()

# Filter genes with standard deviation > 1.5
expr <- expr[standard.deviations > 1.5,]
dim(expr)
# Output will show the dimensions after filtering by standard deviation.

# Produce PCA plot
pca <- prcomp(t(expr), scale. = TRUE)
pca_data <- data.frame(Sample = rownames(pca$x), 
                       PC1 = pca$x[, 1], 
                       PC2 = pca$x[, 2], 
                       Class = classes)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Class)) +
  geom_point() +
  theme_minimal() +
  labs(title = "PCA Plot of RNA-Seq Data", x = "PC1", y = "PC2")

# Question: Based on the PCA plot, do you think it will be easy or hard to classify these samples using machine learning?
# Answer: 
# Based on the PCA plot, we can observe if the samples form distinct clusters based on class labels.
# If the classes overlap significantly in the PCA space, it suggests that classification may be difficult.
# On the other hand, if the classes form distinct clusters, classification is likely to be easier.
