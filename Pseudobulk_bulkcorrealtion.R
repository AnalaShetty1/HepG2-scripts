# Load necessary library
library(dplyr)

# Step 1: Read the CSV files
df1 <- read.csv("xu_4cell.csv", stringsAsFactors = FALSE)
df2 <- read.csv("4cell.csv", stringsAsFactors = FALSE)

# Step 2: Check the structure of the data frames to ensure they have the expected columns
cat("Structure of df1:\n")
str(df1)

cat("\nStructure of df2:\n")
str(df2)

# Step 3: Merge the two data frames on the columns chrom, start, and end
merged_df <- inner_join(df1, df2, by = c("chr", "start", "end"), suffix = c("_file1", "_file2"))

# Step 4: Check the merged data to ensure it looks correct
cat("\nFirst few rows of merged data:\n")
head(merged_df)

# Step 5: Remove rows with NA values in either of the value columns
merged_df_clean <- merged_df %>%
  filter(!is.na(value_file1) & !is.na(value_file2))
if (nrow(merged_df_clean) > 0) {
  cor_result <- cor.test(merged_df_clean$value_file1, merged_df_clean$value_file2, method = "spearman")
  
  # Print the Spearman correlation results
  cat("Spearman correlation: ", cor_result$estimate, "\n")
  cat("P-value: ", cor_result$p.value, "\n")
} else {
  cat("No overlapping data found to compute correlation.\n")
}

# Load the library
library(ggplot2)
# Assuming you've already merged the data and have a clean data frame (merged_df_clean)

# Create a scatter plot with ggplot
ggplot(merged_df_clean, aes(x = value_file1, y = value_file2)) +
  geom_point(color = "blue", alpha = 0.05) +            # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add regression line (no confidence interval shading)
  labs(title = "Scatter Plot of Values with Spearman Correlation",
       x = "Pseudobulk RT", y = "Bulk RT") +
  annotate("text", x = max(merged_df_clean$value_file1), y = max(merged_df_clean$value_file2), 
           label = paste("Spearman Correlation = ", round(cor(merged_df_clean$value_file1, merged_df_clean$value_file2, method = "spearman"), 2)),
           hjust = 1, vjust = 1, color = "black", size = 5) +
  theme_minimal()


                                      