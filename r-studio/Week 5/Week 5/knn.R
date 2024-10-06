install.packages("recipes")
library(recipes)
library(caret)
set.seed(3)

idxTrain <- createDataPartition(y = iris$Species, p = 0.5, list = FALSE)
iris.train <- iris[idxTrain,]
iris.test <- iris[-idxTrain,]