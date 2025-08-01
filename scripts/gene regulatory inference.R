if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("RTN", "org.Hs.eg.db"))
library(RTN)
library(org.Hs.eg.db)
#Prepare Your Data
# Expression matrix: genes in rows, samples in columns
# TPM or counts (normalized)
data_matrix <- as.matrix(df)  # use your UC expression matrix here

# Make sure rownames are gene symbols and unique
rownames(data_matrix) <- make.unique(rownames(data_matrix))

# Filter low-expressed genes (optional)
keep <- rowMeans(data_matrix) > 1
data_matrix <- data_matrix[keep, ]
#Define Transcription Factors
# Download TF list
library(dorothea)
data(dorothea_hs, package = "dorothea")
BiocManager::install("RTN")
library(RTN)

# Filter for high confidence (A and B levels)
tf_regulons <- dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B"))

head(tf_regulons)
tfs <- intersect(unique(tf_regulons$tf), rownames(data_matrix))
#Initialize and Run RTN Analysis
# Create TNI object
tni <- tni.constructor(expData = data_matrix, regulatoryElements = tfs)

# Preprocessing
tni <- tni.preprocess(tni)

# Mutual Information (MI) analysis
tni <- tni.permutation(tni)

# Bootstrap and DPI filter (refine network)
tni <- tni.bootstrap(tni)
tni <- tni.dpi.filter(tni)
# visualization
# Extract the GRN
regulons <- tni.get(tni, what = "regulons.and.mode")
# View one TF’s targets
regulons$TP53  # replace with any TF you're interested in
# Export network as edge list
edge_list <- tni.get(tni, what = "refnet")
write.table(edge_list, "UC_GRN_edge_list.txt", sep="\t", quote=FALSE, row.names=FALSE)
# Visualize (optional)
library(igraph)
g <- graph_from_data_frame(edge_list)
plot(g, vertex.size = 5, vertex.label.cex = 0.6)
tp53_targets <- tni.get(tni, what = "regulons")[["TP53"]]
# Simple visualization using igraph
library(igraph)

# Create edge list: TF → Target
edges <- data.frame(from = "TP53", to = tp53_targets)

# Create graph
g <- graph_from_data_frame(edges, directed = TRUE)

# Plot with target labels
plot(g,
     vertex.label = V(g)$name,
     vertex.size = 20,
     vertex.color = ifelse(V(g)$name == "TP53", "tomato", "skyblue"),
     edge.arrow.size = 0.5,
     layout = layout_with_fr)
# to identify regulators for specific gene (il6 for example due to its centrality)
library(dorothea)
library(tidyverse)

# Load high-confidence TF-target interactions
data(dorothea_hs, package = "dorothea")

# Filter for IL6 targets
il6_regulators <- dorothea_hs %>%
  filter(target == "IL6", confidence %in% c("A", "B"))

# View result
il6_regulators


# df is your TPM expression matrix (genes x samples)
# Transpose your expression matrix so samples are rows and genes are columns
df_t <- t(df)

# Make sure IL6 is in the data
"IL6" %in% colnames(df_t)

# Subset IL6 expression and TF expressions (now as columns)
il6_expr <- df_t[, "IL6"]
tf_exprs <- df_t[, intersect(dorothea_tfs, colnames(df_t))]

# Compute correlation
cor_with_il6 <- apply(tf_exprs, 2, function(tf_expr) cor(tf_expr, il6_expr, method = "pearson"))

# Sort and show top TFs
top_tf_regulators <- sort(cor_with_il6, decreasing = TRUE)
head(top_tf_regulators, 10)

intersection <-  intersect(names(top_tf_regulators),il6_regulators$tf)
names(top_tf_regulators)

# visualize 
library(igraph)
library(ggraph)
library(tidygraph)
top_n <- 10

# Create edge list: TF -> IL6
edges <- data.frame(
  from = names(top_tf_regulators)[1:top_n],
  to = rep("IL6", top_n),
  weight = top_tf_regulators[1:top_n]
)

# Make graph
graph <- tbl_graph(edges = edges, directed = TRUE)

# Plot
ggraph(graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = abs(weight), edge_color = weight > 0),
                 arrow = arrow(length = unit(3, 'mm')),
                 end_cap = circle(4, 'mm')) +
  geom_node_point(size = 6, color = "skyblue") +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_color_manual(values = c("red", "forestgreen"),
                          labels = c("Negative", "Positive"),
                          name = "Correlation") +
  theme_void() +
  ggtitle("Top TFs Regulating IL6 (Correlation-Based)")

# using GENIE3
BiocManager::install("GENIE3")
library(GENIE3)
library(dorothea)
data(dorothea_hs, package = "dorothea")

# Get high-confidence human TFs
tf_list <- unique(dorothea_hs$tf[dorothea_hs$confidence %in% c("A", "B")])

# Filter TFs that exist in your matrix
tfs_in_data <- intersect(tf_list, rownames(df))

# run the model
# Run GENIE3
weight_matrix <- GENIE3(exprMatrix = as.matrix(df), regulators = tfs_in_data, nCores = 4)
link_list <- getLinkList(weight_matrix)
head(link_list)
top_links <- link_list[1:1000, ]  # top 1000 strongest edges
length(link_list$regulatoryGene)
library(igraph)
g <- graph_from_data_frame(top_links, directed = TRUE)
plot(g,
     vertex.label = V(g)$name,  # show gene names
     vertex.label.cex = 0.7,   # label size
     vertex.label.color = "black",
     vertex.size = 5,
     edge.arrow.size = 0.3,
     main = "Top Regulatory Network (GENIE3)")
library(circlize)
# Keep only needed columns
links <- top_links[, c("regulatoryGene", "targetGene", "weight")]

# Optional: Filter top edges if too many, e.g. top 200 by weight
links <- links[order(-links$weight), ][1:200, ]

# Clean links data
links <- na.omit(links)    
# Remove rows with NAs
links$regulatoryGene <- as.character(links$regulatoryGene)
links$targetGene <- as.character(links$targetGene)
links <- links[links$regulatoryGene != links$targetGene, ]  # Remove self-links


# Define colors for genes: regulators tomato, targets skyblue
all_genes <- unique(c(links$regulatoryGene, links$targetGene))
gene_colors <- rep("skyblue", length(all_genes))
names(gene_colors) <- all_genes
gene_colors[links$regulatoryGene] <- "tomato"

# Plot chord diagram
chordDiagram(
  x = links[, 1:2],
  grid.col = gene_colors,
  transparency = 0.5,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  diffHeight  = -0.04,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.05)
)

# Add gene labels around the circle
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
  },
  bg.border = NA
)

# to visualize regulators for IL-6 (recommend than correlation based method)
# Make sure columns are character
links$regulatoryGene <- as.character(links$regulatoryGene)
links$targetGene <- as.character(links$targetGene)

# Filter edges where IL6 is the target
il6_regulators <- link_list[link_list$targetGene == "IL6", ]

# Sort by weight (descending)
il6_regulators <- il6_regulators[order(-il6_regulators$weight), ]

# View top regulators
head(il6_regulators, 10)

# extract the regulators 
library(igraph)

# Optionally, select top N regulators by weight or all (ordered)
top_il6_regulators <- il6_regulators[order(-il6_regulators$weight), ]

# Create an edge list for igraph
edges_il6 <- data.frame(
  from = top_il6_regulators$regulatoryGene,
  to = top_il6_regulators$targetGene
)

# Create graph object
g_il6 <- graph_from_data_frame(edges_il6, directed = TRUE)
plot(g_il6,
     vertex.size = 30,
     vertex.label.color = "black",
     vertex.color = "skyblue",
     edge.arrow.size = 0.5,
     main = "IL6 Regulators Network")
write.csv(edges_il6, "il6_regulators_network.csv", row.names = FALSE)


# Filter IL6 regulators edges
il6_regulators <- link_list %>% filter(targetGene == "IL6")

library(igraph)
library(ggraph)
library(tidyverse)

# Assuming il6_regulators dataframe has columns: regulatoryGene, targetGene, weight
library(igraph)
library(ggraph)
library(tidyverse)

# Assume g is your igraph object with edge attribute 'weight' and vertex attribute 'name'

# Create a color group to highlight IL6
V(g)$color_group <- ifelse(V(g)$name == "IL6", "IL6", "other")

# Plot without edge labels
ggraph(g, layout = "fr") + 
  geom_edge_link(aes(width = weight),    # edge width proportional to weight
                 arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm'),
                 color = "grey50") +
  geom_node_point(aes(color = color_group), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_width(range = c(0.2, 4)) +
  scale_color_manual(values = c("IL6" = "red", "other" = "steelblue")) +
  theme_void() +
  ggtitle("IL6 Regulators Network (no edge labels)")


# Assuming your edge list is called il6_regulators or edges_il6
write.csv(il6_regulators, "IL6_regulators_edges.csv", row.names = FALSE)

# Or tab-delimited (preferred by Cytoscape)
write.table(il6_regulators, "IL6_regulators_edges.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# plot regulators on ridge plot 
# Install if needed
install.packages("ggridges")
library(ggridges)
library(ggplot2)
library(tidyr)
library(dplyr)

# Example regulators list (replace with your actual regulators)
regulators <- il6_regulators$regulatoryGene

# Subset expression matrix to these regulators
expr_regulators <- df[rownames(df) %in% regulators, ]

# Convert to long format: gene, sample, expression
df_long <- expr_regulators %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression")

# Plot ridge plot
ggplot(df_long, aes(x = expression, y = gene, fill = gene)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5) +  # increase spacing here
  theme_minimal() +
  labs(title = "Expression Distribution of Regulators",
       x = "Expression",
       y = "Regulator Gene") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))  # smaller y labels
# smaller text size









