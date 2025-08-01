if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("MuSiC")
devtools::install_github("dviraran/xCell")
library(xCell)

# Make sure your matrix is genes in rows, samples in columns
xcell_res <- xCellAnalysis(as.matrix(df))
head(xcell_res)
# heatmap
library(pheatmap)

pheatmap(xcell_res,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "xCell Cell Type Enrichment Heatmap")
library(reshape2)
library(ggplot2)
library(reshape2)
library(ggplot2)

# Reshape to long format
xcell_long <- melt(xcell_res)
colnames(xcell_long) <- c("CellType", "Sample", "Score")

# bar plot
ggplot(xcell_long, aes(x = Sample, y = Score, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("xCell Cell Type Enrichment by Sample")
# variable only

top_ct <- head(order(apply(xcell_res, 1, var), decreasing = TRUE), 10)
xcell_top <- xcell_res[top_ct, ]

xcell_long <- melt(xcell_top)
colnames(xcell_long) <- c("CellType", "Sample", "Score")
# bar
ggplot(xcell_long, aes(x = Sample, y = Score, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Cell Types Enriched by xCell")

# violin
cell_type <- "aDC"

violin_df <- data.frame(
  Sample = colnames(xcell_res),
  Score = as.numeric(xcell_res[cell_type, ]),
  Condition = coldata$condition
)


# Ensure colnames(xcell_res) match rownames(coldata)
all(colnames(xcell_res) == rownames(coldata))  # Should be TRUE

library(dplyr)

# Build dataframe for violin plot
violin_df <- data.frame(
  Sample = colnames(xcell_res),
  Score = xcell_res[cell_type, ],
  Condition = coldata$condition
)


library(ggplot2)

ggplot(violin_df, aes(x = Condition, y = Score, fill = Condition)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  theme_minimal() +
  ggtitle(paste("xCell Enrichment of", cell_type)) +
  ylab("Enrichment Score") +
  xlab("Condition") +
  theme(plot.title = element_text(hjust = 0.5))

install.packages("ggpubr")
library(ggpubr)

ggplot(violin_df, aes(x = Condition, y = Score, fill = Condition)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Control", "Inactive_UC"),
                                        c("Control", "Active_UC"),
                                        c("Inactive_UC", "Active_UC"))) +
  theme_minimal() +
  ggtitle(paste("xCell Enrichment of", cell_type)) +
  ylab("Enrichment Score") +
  xlab("Condition")

library(ggplot2)

# Specify the cell type you're interested in
cell_type <- "aDC"  # Change to any cell type present in rownames(xcell_res)
# violin and box
# Construct the dataframe
violin_df <- data.frame(
  Sample = colnames(xcell_res),
  Score = as.numeric(xcell_res[cell_type, ]),
  Condition = coldata$condition
)

# Create the violin plot with boxplot inside
ggplot(violin_df, aes(x = Condition, y = Score, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +       # Violin shape
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +  # Boxplot inside
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +        # Optional: add points
  labs(
    title = paste("xCell Score for", cell_type),
    x = "Condition",
    y = "xCell Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.position = "none"
  )

#
ggplot(violin_df, aes(x = Condition, y = Score, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1.2, alpha = 0.5) +
  labs(title = paste("Boxplot of", cell_type, "Score"), y = "Score") +
  theme_minimal()
#
ggplot(violin_df, aes(x = Sample, y = Score, group = Condition, color = Condition)) +
  geom_line() +
  geom_point() +
  labs(title = paste("xCell Score Trend -", cell_type), y = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))




#

library(dplyr)
# Load required package
library(ggplot2)

# Create violin + box plot
p <- ggplot(violin_df, aes(x = Condition, y = Score, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +  # Embedded boxplot
  theme_minimal() +  # Clean theme
  labs(title = paste("Distribution of", cell_type, "Score by Condition"),
       y = "xCell Score",
       x = "Condition") +
  theme(plot.title = element_text(hjust = 0.5))  # Center title

# Print plot (required in scripts/functions)
print(p)

#pheatmap
library(pheatmap)
annotation_df <- data.frame(Condition = coldata$condition)
rownames(annotation_df) <- colnames(xcell_res)  # very important!

# Subset xCell result if needed, or just use full matrix
pheatmap(as.matrix(xcell_res),
         scale = "row",
         annotation_col = annotation_df,
         main = "xCell Scores Heatmap")


# bar

library(dplyr)

bar_df <- violin_df %>%
  group_by(Condition) %>%
  summarise(mean_score = mean(Score), sd_score = sd(Score))

ggplot(bar_df, aes(x = Condition, y = mean_score, fill = Condition)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score), width = 0.2) +
  labs(title = paste("Mean ± SD of", cell_type, "Score"), y = "Mean xCell Score") +
  theme_minimal()


#
library(reshape2)

# Melt the matrix: converts it from wide to long format
long_df <- melt(xcell_res)

# Rename the columns: Var1 = CellType, Var2 = Sample, value = Score
colnames(long_df) <- c("CellType", "Sample", "Score")
# Add condition to the melted data
long_df$Condition <- coldata$condition[match(long_df$Sample, rownames(coldata))]

# Subset a specific cell type
cell_type <- "CD8+ T-cells"  # or any available in your xcell_res
subset_df <- long_df[long_df$CellType == cell_type, ]

# Violin + boxplot
ggplot(subset_df, aes(x = Condition, y = Score, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(title = paste("xCell Score of", cell_type), y = "Score") +
  theme_minimal()
#

library(reshape2)
library(dplyr)
library(ggplot2)

# Melt xcell matrix
long_df <- melt(xcell_res)
colnames(long_df) <- c("CellType", "Sample", "Score")

# Add Condition column by matching Sample to coldata
long_df$Condition <- coldata$condition[match(long_df$Sample, rownames(coldata))]

# Summarize: mean and SD per cell type per condition
dot_df <- long_df %>%
  group_by(CellType, Condition) %>%
  summarise(mean_score = mean(Score), sd_score = sd(Score), .groups = "drop")

#dot plot
ggplot(dot_df, aes(x = Condition, y = CellType)) +
  geom_point(aes(size = mean_score, color = mean_score)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  scale_size(range = c(2, 8)) +
  labs(
    title = "xCell Score Dot Plot",
    x = "Condition",
    y = "Cell Type",
    color = "Mean Score",
    size = "Mean Score"
  ) +
  theme_minimal()









