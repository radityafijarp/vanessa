iris.kmeans <- kmeans(iris.data, 3)
table(iris.kmeans$cluster, iris.species, dnn = c("Cluster", "Class"))
iris.kmeans$centers
iris.data %>% 
  ggplot(aes(x = Sepal.Width, y = Sepal.Length, 
             col = as.factor(iris.kmeans$cluster), shape = iris.species)) +
  scale_color_manual(values = c("darkgreen", "blue", "red")) +
  labs(col = "Cluster", shape = "Species") +
  geom_point() +
  annotate("point", x = iris.kmeans$centers[1,2], y = iris.kmeans$centers[1,1], 
           col="darkgreen", size = 2, shape = 8) +
  annotate("point", x = iris.kmeans$centers[2,2], y = iris.kmeans$centers[2,1], 
           col="blue", size = 2, shape = 8) +
  annotate("point", x = iris.kmeans$centers[3,2], y = iris.kmeans$centers[3,1], 
           col="red", size = 2, shape = 8)

# Number of clusters (k) is determined by the unique stages in aspwood.stages
k <- length(unique(aspwood.stages))

# Perform k-means clustering on the transposed data (samples are rows now)
aspwood.kmeans <- kmeans(t(aspwood), centers = k)

# Check the clustering results
table(aspwood.kmeans$cluster, aspwood.stages, dnn = c("Cluster", "AspWood Stage"))

# Install required packages if not already installed

library(ggplot2)

# Perform PCA on the transposed data to reduce dimensions
aspwood.pca <- prcomp(t(aspwood), scale = TRUE)

# Create a data frame for plotting
pca_data <- data.frame(PC1 = aspwood.pca$x[,1], 
                       PC2 = aspwood.pca$x[,2], 
                       Cluster = as.factor(aspwood.kmeans$cluster),
                       Stage = aspwood.stages)

# Plot the PCA results with clusters
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster, shape = Stage)) +
  geom_point(size = 3) +
  labs(title = "k-means Clustering of AspWood Data",
       x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()

