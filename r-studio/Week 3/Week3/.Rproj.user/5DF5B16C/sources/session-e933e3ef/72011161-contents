# Set the working directory
setwd("D:\\New-folder\\vanessa\\r-studio\\Week 3\\htseq")
#D:\\New-folder\\vanessa\\r-studio\\Week 3\\htseq
# Load necessary libraries
library(reshape2)
library(tools)
library(dplyr)
library(ggplot2)

# Get a list of all .txt files in the folder
file_list <- list.files(pattern = "\\.txt$", full.names = TRUE)

# Initialize an empty dataframe
combined_df <- data.frame()

# Loop through each file
for (file in file_list) {
  # Read the data without headers
  data <- read.table(file, header = FALSE, col.names = c("Gene", "Value"), stringsAsFactors = FALSE)
  
  # Create a sample name from the file name (removing the path and extension)
  sample_name <- file_path_sans_ext(basename(file))
  
  # Add a column for the sample name
  data$Sample <- sample_name
  
  # Combine with the main dataframe
  combined_df <- rbind(combined_df, data)
}

# Optionally, you can view the combined dataframe
head(combined_df)

# Transform the dataframe into a wide format
wide_df <- dcast(combined_df, Gene ~ Sample, value.var = "Value")

# Convert the dataframe to a matrix
gene_matrix <- as.matrix(wide_df[,-1])  # Exclude the first column (Gene names)
rownames(gene_matrix) <- wide_df$Gene  # Set row names as Gene names

# Filter rows where the row names start with "ENS"
filtered_matrix <- gene_matrix[startsWith(rownames(gene_matrix), "ENS"), ]

# Rename the matrix to 'count'
counts <- filtered_matrix

# View the resulting matrix
head(counts)

####################################
library(DESeq2)
# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(condition = factor(c("autism", "autism", "autism", "control", "control", "control"), 
                                          levels = c("control", "autism"))),
  design = ~ condition
)

# Perform differential expression analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Question: How many genes with expression were tested?
num_genes_tested <- nrow(res)
num_genes_tested

# Question: How many genes are DEGs? (padj < 0.05)
res0.05 <- res %>% as.data.frame() %>% filter(padj < 0.05) %>% arrange(padj)
num_DEGs <- nrow(res0.05)
num_DEGs

# Question: How many DEGs are up and how many are down in autism compared to control?
up_DEGs <- nrow(filter(res0.05, log2FoldChange > 0))
down_DEGs <- nrow(filter(res0.05, log2FoldChange < 0))
up_DEGs
down_DEGs

# Question: How many DEGs do we expect to be false positives given the FDR cutoff?
expected_false_positives <- num_DEGs * 0.05
expected_false_positives

# MA plot
plotMA(res, ylim = c(-5, 5))
title("MA Plot")

# Vulcano plot
res %>%
  as.data.frame() %>%
  mutate(log_pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = log2FoldChange, y = log_pvalue)) +
  geom_point(aes(color = padj < 0.05), size = 1) +
  scale_color_manual(values = c("grey", "blue")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(p-value)")

# Create a new table with fold-change cutoff
res0.05FC <- res0.05 %>%
  filter(log2FoldChange > log2(2) | log2FoldChange < log2(0.5))

# Save the results for later
save(dds, vst, res, res0.05, res0.05FC, file = "autism_DEA.RData")

# GO enrichment analysis
# Convert ensembl gene IDs to gene symbols for enrichr
conversion_table <- read_tsv("Human_ensembl_ids_to_symbols.txt")
DEG_symbols <- res0.05 %>%
  rownames_to_column("ensembl") %>%
  left_join(conversion_table, by = "ensembl") %>%
  select(symbol) %>%
  filter(!is.na(symbol))

# Write gene symbols to file for enrichr input
write_lines(DEG_symbols$symbol, "DEG_symbols.txt")

