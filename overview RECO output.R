# Import necessary packages
library(dplyr)
library(ggplot2)
library(readr)
library(reshape2)
library(viridis)
library(tidyr)
library(stringr)
library(RColorBrewer)

# List all the target_sequence_counts files
files <- list.files(pattern = "*_counts.csv")

# Sort files to get the first file and remaining files
files <- sort(files)
first_file <- files[1]
remaining_files <- files[-1]  # Excludes the first element

# Read in the first file
df <- read_csv(first_file)

# Read in the remaining files (without headers)
for(file in remaining_files){
  temp_df <- read_csv(file, skip = 1, col_names = FALSE) 
  # Make sure the column names match
  colnames(temp_df) <- colnames(df)
  
  df <- bind_rows(df, temp_df)
}

# Save the concatenated DataFrame to a new csv file
write_csv(df, "combined_counts.csv")

## Read the data
df2 <- read_csv("combined_counts.csv")

# Remove the "_multi_variants.fasta" part from the Species and Strain columns
df2$Filename <- str_replace(df2$Filename, "_multi_variants.fasta", "")

# Melt the DataFrame
df_melt <- melt(df2, id.vars = c("Filename", "Target Sequence"))


# Plot
# Create the bar plot
ggplot(df_melt, aes(x = `Target Sequence`, y = value)) +
  geom_point(aes(color = Filename)) +
#  scale_color_brewer(palette = "Dark2") +
  theme(strip.text.x = element_text(angle = 0, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "bottom") +
  scale_y_log10() +
  scale_color_viridis_d(direction = 1)+
  facet_grid( ~ variable)
                                   