library(readr)
#install.packages("tibble")  # Install if you haven't already
library(tibble)
library(dplyr)

# GO enrichment analysis
setwd("D:\\New-folder\\vanessa\\r-studio\\Week 3\\Week3")
# Convert ensembl gene IDs to gene symbols for enrichr
conversion_table <- read_tsv("Human_ensembl_ids_to_symbols.txt")
DEG_symbols <- res0.05 %>%
  rownames_to_column("ensembl_id") %>%
  left_join(conversion_table, by = "ensembl_id") %>%
  select(gene_symbol) %>%
  filter(!is.na(gene_symbol))

# Write gene symbols to file for enrichr input
write_lines(DEG_symbols$gene_symbol, "DEG_symbols.txt")
