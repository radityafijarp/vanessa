if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install the DESeq2 package
BiocManager::install("DESeq2")

library(tidyverse)

read_tsv("AspWood/K1-01_sortmerna_trimmomatic/quant.sf") %>% 
  slice(1:10)

samples <- list.files("AspWood")

counts <- tibble()
for (sample in samples) {
  file <- paste0("AspWood/", sample, "/quant.sf")
  
  sample.trimmed <- gsub("_sortmerna_trimmomatic", "", sample)
  
  c <- read_tsv(file) %>%
    select(Name, NumReads) %>%
    rename(Genes = Name, !!sym(sample.trimmed) := NumReads)
  
  if (sample == samples[1]) {
    counts <- c 
  } else {
    counts <- cbind(counts, c %>% select(-Genes))
  }
}

dim(counts)

counts %>% 
  select(Genes, `K1-01`) %>% 
  filter(Genes == "Potra2n765s36714.1") 

# Method 1
counts <- counts %>% 
  pivot_longer(cols = -Genes, names_to = "Samples", values_to = "Expression") %>%
  separate(Genes, into = c("Genes"), sep = "\\.", extra = "drop") %>% 
  group_by(Genes, Samples) %>%
  summarise(Expression = sum(Expression)) %>%
  pivot_wider(names_from = "Samples", values_from = "Expression")

# Method 2
counts <- counts %>% 
  separate(Genes, into = c("Genes"), sep = "\\.", extra = "drop") %>% 
  group_by(Genes) %>%
  summarise_if(is.numeric, sum)

dim(counts)

#
counts.long <- counts %>% 
  pivot_longer(cols = -Genes, names_to = "Samples", values_to = "Expression")

dim(counts.long)

#Boxplot of Counts
library(ggplot2)

# Boxplot without outliers
counts.long %>%
  ggplot(aes(x = Samples, y = Expression)) +
  geom_boxplot(outliers = FALSE) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Create a new boxplot where the expression values have been log2-transformed. Use log2(counts+1)! Why?
# Log-transformed boxplot
counts.long %>%
  mutate(Expression = log2(Expression + 1)) %>%
  ggplot(aes(x = Samples, y = Expression)) +
  geom_boxplot(outliers = FALSE) 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Effect of Log-transformation:
counts.long %>%
  ggplot(aes(x = Expression)) +
  geom_histogram(bins = 50) +
  ggtitle("Histogram of Counts")


counts.long %>%
  mutate(logExpression = log2(Expression + 1)) %>%
  ggplot(aes(x = logExpression)) +
  geom_histogram(bins = 50) +
  ggtitle("Histogram of Log-transformed Counts")

##Normalization
#library(DESeq2)

counts.mat <- counts %>% 
  column_to_rownames(var = "Genes") %>% 
  as.matrix() %>% 
  round()

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts.mat,
                                      colData = data.frame(conditions = as.factor(colnames(counts.mat))),
                                      design = ~ conditions)

#Size Factor Normalization:
dds <- DESeq2::estimateSizeFactors(dds)

DESeq2::counts(dds, normalized=TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Genes") %>%
  pivot_longer(cols = -Genes, names_to = "Samples", values_to = "Expression") %>%
  mutate(Expression = log2(Expression + 1)) %>%
  ggplot(aes(x = Samples, y = Expression)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Comparison of Pre- and Post-Normalization Boxplots:
size_factors <- DESeq2::sizeFactors(dds)
print(size_factors)

###Variance Stabilizing Transformation
# Applying variance stabilizing transformation
vst <- DESeq2::varianceStabilizingTransformation(counts.mat)

# Subtract the minimum value to ensure that genes without mapped reads have value zero
vst <- vst - min(vst)

#Create a boxplot of the VST-normalized values for each sample, as before
# Convert VST-normalized data into long format for ggplot
vst.long <- as.data.frame(vst) %>%
  rownames_to_column(var = "Genes") %>%
  pivot_longer(cols = -Genes, names_to = "Samples", values_to = "Expression")

# Create boxplot
vst.long %>%
  ggplot(aes(x = Samples, y = Expression)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

##Removing Lowly Expressed Genes:
#Genes not expressed at all:
#How many genes are not expressed at all in our data?
# Count the number of genes with zero expression across all samples
no_expression <- sum(rowSums(vst == 0) == ncol(vst))
no_expression

#How many genes have less than 10 mapped reads in total in our data? Remove these genes from vst
#Genes with fewer than 10 mapped reads in total:
# Count the number of genes with less than 10 reads across all samples
low_expression <- sum(rowSums(vst) < 10)
low_expression

# Filter out genes with fewer than 10 mapped reads
vst_filtered <- vst[rowSums(vst) >= 10, ]

#Save the data for later
save(dds, vst, file="AspWood_normalized.RData")

# Plot expression of SUS6 (Potra2n4c9149)
vst.long %>%
  filter(Genes == "Potra2n4c9149") %>%
  ggplot(aes(x = Samples, y = Expression)) +
  geom_line() +
  geom_point() +
  ggtitle("Expression Profile of SUS6 (Potra2n4c9149)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plot expression of SUS6 (Potra2n4c9149) - Line chart without points
vst.long %>%
  filter(Genes == "Potra2n16c30563") %>%
  ggplot(aes(x = Samples, y = Expression, group = 1)) +  # 'group = 1' ensures the data is treated as a continuous line
  geom_line() +
  ggtitle("Expression Profile of SUS6 (Potra2n16c30563)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#Plot expression of SUS6 (Potra2c131s34676) - Line chart without points
vst.long %>%
  filter(Genes == "Potra2c131s34676") %>%
  ggplot(aes(x = Samples, y = Expression, group = 1)) +  # 'group = 1' ensures the data is treated as a continuous line
  geom_line() +
  ggtitle("Expression Profile of SUS6 (Potra2c131s34676)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#Run all gene

# Get unique genes
unique_genes <- unique(vst.long$Genes)

# Loop over each gene
for (i in unique_genes) {
  # Create the plot for the current gene
  p <- vst.long %>%
    filter(Genes == i) %>%
    ggplot(aes(x = Samples, y = Expression, group = 1)) +  # 'group = 1' ensures the data is treated as a continuous line
    geom_line() +
    ggtitle(paste("Expression Profile of", i)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  # Save the plot to a file (e.g., PNG format)
  ggsave(filename = paste0("Expression_Profile_", i, ".png"), plot = p, width = 8, height = 6)
}
