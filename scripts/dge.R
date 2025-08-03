library(here)
here()
install.packages("renv")
renv::init()
load("project_files/data/matrix.RData")
load("project_files/data/data.RData")
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
library(gt)
library(gtsummary)
library(chromote)

mytable <- metadata %>%
  mutate(three_year_survival= factor(three_year_survival, levels= c(0,1), labels = c("alive", "dead"))) %>%
  select(race, gender, ethnicity, prior_malignancy, year_of_diagnosis, age_at_diagnosis, survival_days, three_year_survival ) %>%
  tbl_summary(
    by = three_year_survival, missing_text =  "missing"
  ) %>%
  as_gt()
gt::gtsave(mytable, filename = "characteristics.png")
library(caret)
library(survival)
install.packages("survminer")
library(survminer)
metadata <- downSample(metadata, metadata$three_year_survival, list = FALSE, yname = "three_year_survival")

surv_obj <- Surv(time = metadata$survival_days, event = metadata$survival_status)

fit <- survfit(surv_obj ~ 1)  # Overall survival
ggsurvplot(fit,
           data = metadata,
           conf.int = TRUE,
           risk.table = TRUE,
           title = "Overall Survival",
           xlab = "Days",
           ylab = "Survival Probability")
library(ggplot2)

install.packages("ggpubr")
library(ggpubr)

# Plot
table(metadata$three_year_survival, metadata$ethnicity)
chi_res <- chisq.test(table(metadata$three_year_survival, metadata$gender))
p_val <- chi_res$p.value

metadata %>%
  mutate(three_year_survival= factor(three_year_survival, levels= c(0,1), labels = c("alive", "dead"))) %>%
  ggplot(aes(x = three_year_survival, fill = gender)) +
  geom_bar(position = "dodge") +
  labs(
    title = "3-Year Survival",
       x = "survival status",
       y = "gender", 
      fill = "gender"
  ) +
  theme_minimal() +
  annotate("text", x = 1.5, y = max(table(metadata$three_year_survival)), 
           label = paste0("Chi-square p = ", signif(p_val, 3)), 
           size = 5, fontface = "italic")

#


library(ggplot2)
metadata %>%
  mutate(three_year_survival= factor(three_year_survival, levels= c(0,1), labels = c("alive", "dead"))) %>%  
  ggplot(aes(x = three_year_survival, y = age_at_diagnosis, fill = three_year_survival)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Survival Days by 3-Year Outcome",
       x = "3-Year Survival",
       y = "Survival Days") +
  theme_minimal() +
  theme(legend.position = "lowerbottom")

metadata %>%
  mutate(survival_status = factor(survival_status, levels = c(0, 1), labels = c("alive", "dead"))) %>%  
  ggviolin(x = "gender", y = "survival_days",
           fill = "survival_status", palette = "jco",
           add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(aes(group = survival_status), method = "t.test") +
  labs(
    title = "Survival Days by Gender and Overall Survival Status",
    x = "Gender",
    y = "Survival Days",
    fill = "Survival Status"
  ) +
  theme_minimal()









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

#
all_genes <- row.names(dds)
all_genes <- unlist(all_genes)
all_genes <- str_replace(all_genes, "\\.\\d+$", "")

# Retrieve mappings
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

all_mapping <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  values = all_genes,
  mart = mart
)
sig_mapping <-  getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  values = sig_res$gene,
  mart = mart
)
# Merge with results
# GO enrihcment 
ego <- enrichGO(gene = sig_mapping$entrezgene_id,
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
kk <- enrichKEGG(gene         = sig_mapping$entrezgene_id,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
library(DOSE)
library(clusterProfiler)
x <- enrichDO(gene          = sig_mapping$entrezgene_id,
              ont           = "HDO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = as.character(all_mapping$entrezgene_id),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(x)
barplot(x)

genelist <- res$log2FoldChange
genelist <- unique(genelist)
names(genelist) <- all_mapping$entrezgene_id

genelist <- sort(genelist, decreasing = TRUE)
genelist <- genelist[!is.na(names(genelist)) & !duplicated(names(genelist))]
ego3 <- gseGO(geneList     = genelist,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
goplot(ego)
#
BiocManager::install("msigdbr")
library(msigdbr)
genesets <- msigdbr()
# Show all pathway names
sets <- msigdbr(species = "Homo sapiens", category = "C2")
unique(sets$gs_collection_name)
# Get KEGG Medicus gene sets for Homo sapiens
Reactome_Pathways <- msigdbr(species = "Homo sapiens")
Reactome_Pathways <- subset(Reactome_Pathways, gs_collection_name == "Reactome Pathways")
term2gene <- Reactome_Pathways[, c("gs_name", "gene_symbol")]
# Unique gene symbols from your term2gene table
gene_symbols <- unique(term2gene$gene_symbol)

gene_list <- res$log2FoldChange
res <- results(dds)

# Extract Ensembl IDs from rownames
ensembl_ids <- rownames(res) |> str_replace("\\.\\d+$", "")

# Create res_df and attach Ensembl IDs
res_df <- as.data.frame(res) %>%
  mutate(ensembl_id = ensembl_ids)

# Join DESeq2 results with Entrez IDs
res_clean <- res_df %>%
  inner_join(all_mapping, by = c("ensembl_id" = "ensembl_gene_id"))


names(res)

  
res_clean <- res_clean[!duplicated(res_clean$ensembl_id),]

row.names(res_clean) <- res_clean$ensembl_id
genes <- res_clean$log2FoldChange
names(genes) <- row.names(res_clean)

# Map to ENTREZ IDs
symbol_to_entrez <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters    = "hgnc_symbol",
  values     = gene_symbols,
  mart       = mart
)

# Remove entries without ENTREZ ID
symbol_to_entrez <- symbol_to_entrez[!is.na(symbol_to_entrez$entrezgene_id), ]
symbol_to_entrez <- symbol_to_entrez[!duplicated(symbol_to_entrez$hgnc_symbol), ]

term2gene_entrez <- term2gene %>%
  inner_join(symbol_to_entrez, by = c("gene_symbol" = "hgnc_symbol")) %>%
  dplyr::select(gs_name, entrezgene_id)
genelist <- sort(genes, decreasing = TRUE)
names(genelist) <- res_clean$entrezgene_id
genelist <- genelist[!is.na(names(genelist)) & !duplicated(names(genelist))]
gsea_result <- GSEA(
  geneList    = genelist,   # named numeric vector: names = SYMBOLs, values = logFC/stat
  TERM2GENE   = term2gene_entrez,
  pvalueCutoff = 0.01
)
length(intersect(names(genelist), term2gene_entrez$entrezgene_id))
enrichplot::dotplot(gsea_result, showCategory = 20)
class(gsea_result)

# This works correctly:
enrichplot::dotplot(gsea_result, showCategory = 20) +
  theme(
    axis.text.y = element_text(size = 5),  # Y-axis (pathway) label size
    axis.text.x = element_text(size = 5),  # X-axis label size
    axis.title = element_text(size = 5),   # X/Y axis title size
    plot.title = element_text(size = 5, face = "bold")  # Optional
  )
