library(ggplot2)
library(ggrepel)  # Load the ggrepel library

# read the distance matrix from the file
dist_matrix <- read.csv("distance_matrix.csv", row.names = 1)

# Convert the DataFrame to a matrix again
dist_matrix <- as.matrix(dist_matrix)

# Convert the matrix to a dist object
dist_obj <- as.dist(dist_matrix)

# Perform MDS
mds_result <- cmdscale(dist_obj)

# Convert the result to a data frame for plotting
mds_df <- data.frame(Label = rownames(mds_result), Dim1 = mds_result[,1], Dim2 = mds_result[,2])

# Create a vector of your unique labels
unique_labels <- paste0("Unique_", 0:61)

# Create a column for the color, set to "red" for unique labels and "black" for all others
mds_df$color <- ifelse(mds_df$Label %in% unique_labels, "red", "lightgrey")

# Create a column for the labels, which should only have values for the unique labels
# Use the sub function to remove "Unique_" from the labels
mds_df$plot_label <- ifelse(mds_df$Label %in% unique_labels, sub("Unique_", "", mds_df$Label), "")

ggplot(data = mds_df, aes(x = Dim1, y = Dim2, color = color, label = plot_label)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_text_repel(size = 2) +  # Adjust size here. Smaller values will yield smaller font size.
  scale_color_identity() +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  ggtitle("MDS of Distance Matrix") +
  theme_minimal()

