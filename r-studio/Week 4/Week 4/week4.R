load("AspWood_normalized.RData")


# Load the libraries
library(tidyverse)
library(usedist)
library(dendextend)
library(pheatmap)

vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>% 
  mutate(Expressed = Expression > 5,
         Silent = Expression < 1) %>%
  group_by(Gene) %>% 
  summarise(Expressed = sum(Expressed),
            Silent = sum(Silent)) %>% 
  filter(Expressed >= 1, Silent >= 1) %>% 
  pull(Gene) -> regulated.genes

aspwood <- vst[rownames(vst) %in% regulated.genes,]

dim(aspwood)

# Plot the expression profile of a random gene
set.seed(3)

aspwood %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Gene") %>% 
  filter(Gene == sample(regulated.genes, 1)) %>% 
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>% 
  separate(Sample, into = c("Tree", "Sample"), sep = "-") %>% 
  mutate(Sample = as.numeric(Sample)) %>% 
  ggplot(aes(x = Sample, y = Expression)) +
  geom_line(linewidth = 1.5) +
  theme_bw()

remove(dds, vst)

#Briefly explain what the code above is doing.
# Read assignment 2, visualize random gene, remove data assignment 2


