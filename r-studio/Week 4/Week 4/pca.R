# Perform PCA on the AspWood data (transpose the data to have samples as rows)
aspwood.pca <- prcomp(t(aspwood), scale. = TRUE)

# View the summary to check the variance explained by each component
summary(aspwood.pca)

# Create a data frame for the PCA results
pca_data <- data.frame(PC1 = aspwood.pca$x[, 1], 
                       PC2 = aspwood.pca$x[, 2], 
                       Stage = aspwood.stages)

# Plot the PCA results
ggplot(pca_data, aes(x = PC1, y = PC2, color = Stage)) +
  geom_point(size = 3) +
  labs(title = "PCA of AspWood Gene Expression Data",
       x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

summary(aspwood.pca)
#Interpretation of Variance
# PC1 explains 43.29% of the variance.
# PC2 explains 22.66% of the variance.
# Together, the first two components (PC1 and PC2) explain 65.95% of the total variance in the data.
# Analysis
# The first two principal components capture a significant portion of the variance, making them suitable for visualizing the major trends in the data. The cumulative proportion tells us that these two components represent a majority of the data's structure, although including more components would further explain more of the variance.
