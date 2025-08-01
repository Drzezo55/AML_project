library(here)
here()
install.packages("renv")
renv::init()
load("data/matrix.RData")
load("data/data.RData")
exp <- expr_df
metadata <- data
View(exp)
BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("apeglm", force = TRUE)
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), force =  TRUE)
BiocManager::install("DOSE", force = TRUE)
install.packages("vsn")
BiocManager::install("MuSiC")
install.packages("remotes")
remotes::install_github("dviraran/xCell")
install.packages("pheatmap")

library(xCell)
library(pheatmap)
library(vsn)
library(DOSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(apeglm)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(stringr)
library(dplyr)
library(caret)

# the ensemble ids has version at its end which will prevent its conversion
ens_id <- str_replace(row.names(exp), "\\.\\d+$", "")
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ens_id,      # ENSG IDs
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

exp_filtered <- exp[!duplicated(exp)& rowSums(is.na(exp)==0),]
exp_filtered$ens_ids <- row.names(exp_filtered)
row.names(exp_filtered) <- NULL
exp_filtered <- exp_filtered[!duplicated(exp_filtered)& rowSums(is.na(exp_filtered)==0),]
row.names(exp_filtered) <- exp_filtered$ens_ids
exp_filtered$ens_ids <- NULL
anyDuplicated(ens_ids)
View(metadata)
metadata <- as.data.frame(colData(metadata))
metadata <- metadata %>%
  mutate(
    survival_status = ifelse(vital_status == "Dead", 1, 0),
    survival_days = days_to_last_follow_up
  )
metadata <- metadata %>%
  mutate (three_year_survival = ifelse(survival_status == 1 & survival_days < 3*365, 1, 0))
metadata$three_year_survival
metadata <- metadata[!is.na(metadata$three_year_survival), ]
summary(metadata$three_year_survival)
metadata$three_year_survival <- as.factor(metadata$three_year_survival)
metadata <- downSample(metadata, metadata$three_year_survival, list = FALSE, yname = "three_year_survival")
metadata$three_year_survival <- NULL

metadata <- metadata %>%
    dplyr::select(barcode, race, gender, age_at_diagnosis, ethnicity,
    survival_status, survival_days, three_year_survival)
row.names(exp_filtered)

df <- exp_filtered[,colnames(exp_filtered) %in% metadata$barcode]
View(df)
#
sum(duplicated(row.names(df)))
dds <- DESeqDataSetFromMatrix(countData=df, 
                              colData=metadata, 
                              design=~three_year_survival)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
res <- results(dds)
dds <- estimateSizeFactors(dds)
res
dds$three_year_survival <- factor(dds$three_year_survival, levels = c("0", "1"), labels = c("alive", "dead"))
# Re-run DESeq with new labels
dds <- DESeq(dds)
# Then extract results
res <- results(dds, contrast = c("three_year_survival", "dead", "alive"))
resultsNames(dds)
summary(res)
res05 <- results(dds, alpha=0.05)
summary(res05)
png("results/MAplot.png", width = 6, height = 5, units = "in", res = 300)
plotMA(res, ylim=c(-2,2))
dev.off()
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="three_year_survival_dead_vs_alive", type="apeglm")

resLFC
resOrdered <- res[order(res$pvalue),]
png("results/MAplot_lfc.png", width = 6, height = 5, units = "in", res = 300)
plotMA(resLFC, ylim=c(-2,2))
dev.off()
# Convert result to dataframe
res_df <- as.data.frame(resLFC)
res_df$gene <- rownames(res_df)

# Filter significant genes
sig_res <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]

# Define up and downregulated
upregulated <- sig_res[sig_res$log2FoldChange > 1, ]
downregulated <- sig_res[sig_res$log2FoldChange < -1, ]

# Save to files
write.csv(upregulated, "results/upregulated_genes.csv", row.names = FALSE)
write.csv(downregulated, "results/downregulated_genes.csv", row.names = FALSE)
res_df$significance <- with(res_df, ifelse(padj < 0.05 & abs(log2FoldChange) > 1,
                                           ifelse(log2FoldChange > 1, "Up", "Down"),
                                           "Not Sig"))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme(legend.title = element_blank())
ggsave("results/volcano_plot.png", width = 6, height = 5, dpi = 300)
#
# If your genes are Ensembl IDs:
sig_res$gene <- row.names(sig_res)
sig_res$gene <- str_replace(sig_res$gene, "\\.\\d+$", "")
genes <- sig_res$gene
gene_mapping <- bitr(genes, fromType = "ENSEMBL",
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

# Merge with results
sig_mapped <- merge(sig_res, gene_mapping, by.x = "gene", by.y = "ENSEMBL")
# GO enrihcment 
ego <- enrichGO(gene = sig_mapped$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",   # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

head(ego)
png("results/go.png", width = 6, height = 5, units = "in", res = 300)
barplot(ego, showCategory = 20)
dev.off()
gene_list <- sig_mapped$log2FoldChange
names(gene_list) <- sig_mapped$ENTREZID
# Sort in decreasing order


