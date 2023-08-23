library(ggplot2)
library(gridExtra)

# File name generator based on perplexity and learning rate
generate_filename <- function(perplexity, learning_rate) {
  return(paste0('tsne_results_CN_20_perplexity_', perplexity, '_lr_', learning_rate, '.csv'))
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
    
    # Generate plot
    p <- ggplot(data, aes(x = Dim1, y = Dim2)) +
      geom_point(aes(color = Strain), size = 2) + 
      labs(title = paste("Perplexity:", perplexity, "LR:", learning_rate)) +
      theme_light() +
      theme(legend.position = "none")
    
    # Append plot to list
    plots[[length(plots) + 1]] <- p
  }
}

# Arrange plots in a grid
grid.arrange(grobs = plots, ncol = 4)
