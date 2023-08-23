library(ggplot2)
library(gridExtra)
library(ggrepel)
library(readr)
library(ggpubr)

setwd("C:/Users/nvanreet/Desktop/Scripts/TSNE")

meta_data <- read.csv("strain_subspecies_ab_meta.csv", sep = ";")
meta_data$Strain <- trimws(meta_data$Strain)

# File name generator based on perplexity and learning rate
generate_filename <- function(perplexity, learning_rate) {
  return(paste0('tsne_results_CN_100_perplexity_', perplexity, '_lr_', learning_rate, '.csv'))
}

# List of perplexities and learning rates
perplexities <- c(10, 30, 50, 70)
learning_rates <- c(100, 200, 500, 1000)

# Generate plots and store them in a list
plots <- list()

for (perplexity in perplexities) {
  for (learning_rate in learning_rates) {
    # Read data
    file_name <- generate_filename(perplexity, learning_rate)
    data <- read.csv(file_name)
    merged_data <- merge(data, meta_data, by = "Strain")
    
    # Generate plot
    p <- ggplot(merged_data, aes(x = Dim1, y = Dim2)) +
      geom_point(aes(color = Abbreviation), size = 2) + 
      labs(title = paste("Perplexity:", perplexity, "LR:", learning_rate)) +
      theme_light() 
    #theme(legend.position = "none")
    
    # Append plot to list
    plots[[length(plots) + 1]] <- p
  }
}

# Arrange plots in a grid with a common legend at the bottom
ggarrange(plotlist = plots, ncol = 4, nrow = 4, common.legend = TRUE, legend = "bottom")

### Select Clusters

tsne <- read_csv("tsne_results_CN_100_perplexity_10_lr_100.csv")
ggplot(tsne, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = Strain), size = 2) + 
  geom_text_repel(aes(label=Strain), size=3) +
  theme_light() +
  theme(legend.position = "none")

# Elbow Method to determine the optimal number of clusters
set.seed(42)
max_k <- 10  # Maximum number of clusters to check
wcss <- numeric(max_k)

for (k in 1:max_k) {
  kmeans_result <- kmeans(tsne[,1:2], centers=k)
  wcss[k] <- kmeans_result$tot.withinss
}

# Plot WCSS to find the 'elbow'
plot(1:max_k, wcss, type="b", pch=19, frame=FALSE, 
     xlab="Number of clusters K", ylab="Total within-clusters sum of squares",
     main="Elbow Method for Determining Optimal K")
lines(lowess(1:max_k, wcss), col="red", lty=2)



# Using kmeans clustering
set.seed(42)  # for reproducibility
number_of_clusters <- 4  # Replace with the number of clusters you expect or determined
km <- kmeans(tsne[,1:2], centers = number_of_clusters)

# Adding the clustering results to the dataframe
tsne$cluster <- as.factor(km$cluster)
merged_data <- merge(tsne, meta_data, by = "Strain")

# Visualizing with ggplot2
library(ggplot2)
ggplot(merged_data, aes(x=Dim1, y=Dim2, shape=cluster, color = Abbreviation)) + 
  geom_point(alpha=0.8, size=3) + 
  geom_text_repel(aes(label=Strain), size=2) +
  theme_minimal() +
  labs(title="t-SNE visualization with K-means Clustering")
