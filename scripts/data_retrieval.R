# Install TCGAbiolinks
BiocManager::install("TCGAbiolinks")
BiocManager::install(c("SummarizedExperiment", "DESeq2", "edgeR",
"ggplot2", "ComplexHeatmap", "clusterProfiler","org.Hs.eg.db"), force = TRUE)
install.packages("DT")
install.packages("groupdata2")
install.packages("caret")      # if not installed
library(caret)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DT)
library(here)
library(ggplot2)
library(dplyr)

project_name <- "TCGA-LAML"
query <- GDCquery(project = project_name,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts")

metadata <- query[[1]][[1]]
GDCdownload(query, method = "api")
data = GDCprepare(query)
save(data, file = "data.RData")

data  # check the structure
data   # expression matrix
colData(data)[1:5, ]   # clinical metadata
class(data)
colData(data)[1:5, ]  # View first 5 samples' metadata
View(colData(data))  # Opens a spreadsheet-like viewer in RStudio
datatable(as.data.frame(colData(data)))
data@rowRanges$gene_id
# View raw gene expression counts
head(assay(data))  # Shows first 6 genes across all samples
expr_df <- as.data.frame(assay(data))
# Gene identifiers (usually Ensembl or gene symbols)
head(rownames(data))
# Sample identifiers (barcodes)
head(colnames(data))

# Example: Boxplot of TP53 expression across sample types
colData(data)
metadata
clinical <- colData(data)
df <- as.data.frame(clinical)
View(df)

df <- df %>%
  mutate(
    survival_status = ifelse(vital_status == "Dead", 1, 0),
    survival_days = days_to_last_follow_up
  )
df <- df %>%
  mutate (three_year_survival = ifelse(survival_status == 1 & survival_days < 3*365, 1, 0))
df$three_year_survival
factors <- c()
summary(df$three_year_survival)
df$three_year_survival <- as.factor(df$three_year_survival)
df$three_year_survival <- as.character(df$three_year_survival)

new_df <- df[!is.na(df$three_year_survival), ]
new_df$three_year_survival <- as.factor(new_df$three_year_survival)
new_df <- downSample(new_df, new_df$three_year_survival, list = FALSE, yname = "three_year_survival")
new_df$three_year_survival
levels(new_df$three_year_survival)
(new_df$three_year_survival)
new_df <- new_df[, !duplicated(names(new_df))]

ggplot(new_df, aes(x = three_year_survival)) +
  geom_bar(fill = 'red') +
  labs(x = '3-Year Survival Status', y = 'Count', title = '3-Year Survival Distribution') +
  theme_minimal()
 new_df <- new_df %>%
  dplyr::select(c(barcode, race, gender, ethnicity, survival_status, 
    survival_days, three_year_survival))

new_df$barcode
summary(new_df)


