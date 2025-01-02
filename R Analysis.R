library(readr)
library(affy)
library(ggfortify)
library(rgl)
library(plot3D)
library(plotly)
library(stats)
library(scatterplot3d)
library(genefilter)
library(matrixStats)
library(ComplexHeatmap)
library(affyPLM)
library(pd.hg.u133.plus.2)
library(readxl)
library(hgu133a.db)
library(hgu133plus2.db)
library(hgu133acdf)
library(hgu133plus2cdf)
library(multtest)
library(EnhancedVolcano)

#
############# import CEL files ###############
celfilesGSE7476_directory = "C:/Users/احمد/OneDrive - Nile University/Desktop/omics project/GSE19804_RAW"
cels_GSE7476 = list.files(celfilesGSE7476_directory, pattern = "CEL")

affydata_GSE7476 = ReadAffy(celfile.path =celfilesGSE7476_directory )
affydata_GSE7476
class(affydata_GSE7476)

#
#metadata = read_table("phenotype data (library(readxl)
#metadata = read_table("phenotype data (library(readxl)
metadata <- read_excel("C:/Users/احمد/OneDrive - Nile University/Desktop/omics project/metadata.xlsx")
View(metadata)

##################### data exploration ############3
class(affydata_GSE7476)
sampleNames(affydata_GSE7476)
featureNames(affydata_GSE7476)
head(featureNames(affydata_GSE7476))
annotation(affydata_GSE7476)
dim(affydata_GSE7476)

# Raw expression matrix
#raw_expression <- exprs(affy_data)
#boxplot(affy_data, main = "Box Plot (Raw Data)", col = rainbow(ncol(affy_data)))
#hist(affy_data, main = "Histogram (Raw Data)")

# Import CEL files
cel_files_directory <- "C:/Users/احمد/OneDrive - Nile University/Desktop/omics project/GSE19804_RAW"
cel_files <- list.files(cel_files_directory, pattern = "CEL")

# Load Affymetrix CEL files
affy_data <- ReadAffy(celfile.path = cel_files_directory)

# Check if the AffyBatch object is loaded
print(affy_data)
class(affy_data)

# Extract raw expression matrix
raw_expression <- exprs(affy_data)
head(raw_expression)

# PCA on raw expression data
raw_expression_t <- t(raw_expression)
pca_raw <- prcomp(raw_expression_t, center = TRUE, scale. = TRUE)
autoplot(pca_raw, data = metadata, colour = 'Condition')

# Preprocessing (Normalization and summarization)
processed_eset <- threestep(affy_data, 
                            background.method = "IdealMM", 
                            normalize.method = "quantile", 
                            summary.method = "average.log")
processed_expression <- exprs(processed_eset)
boxplot(processed_eset, main = "Box Plot (Processed Data)", col = rainbow(ncol(processed_eset)))
hist(processed_eset, main = "Histogram (Processed Data)")

# Save processed data
write.exprs(processed_eset, file = "processed_expression_data.txt")
processed_data <- read_delim("processed_expression_data.txt", delim = "\t")
names(processed_data)[1] <- "Probe_ID"

# Map probe IDs to gene symbols
probe_to_symbol <- as.data.frame(hgu133aSYMBOL)
colnames(probe_to_symbol) <- c("Probe_ID", "Gene_Symbol")
mapped_data <- merge(processed_data, probe_to_symbol, by = "Probe_ID", all.x = TRUE)

# Handle NAs in gene symbols
mapped_data[is.na(mapped_data$Gene_Symbol), "Gene_Symbol"] <- paste0("UNKNOWN_", 1:sum(is.na(mapped_data$Gene_Symbol)))
non_na_data <- mapped_data[!is.na(mapped_data$Gene_Symbol), ]
duplicates <- duplicated(non_na_data$Gene_Symbol)
sum_duplicates <- sum(duplicates)

# Aggregate expression values by gene symbol
aggregated_data <- aggregate(. ~ Gene_Symbol, data = non_na_data[-1], FUN = mean)
rownames(aggregated_data) <- aggregated_data$Gene_Symbol
aggregated_data <- aggregated_data[ , -1]

# Save the final gene expression matrix
save(aggregated_data, file = "final_gene_expression_matrix.RData")
write.csv(aggregated_data, file = "final_gene_expression_matrix.csv")

# PCA on processed data
gene_expression_matrix <- as.matrix(aggregated_data)
pca_processed <- prcomp(t(gene_expression_matrix), center = TRUE, scale. = TRUE)
autoplot(pca_processed, data = metadata, colour = 'Condition')

# Differential Expression Analysis (DEA)

groups <- unique(metadata$Condition)
group1 <- groups[1]
group2 <- groups[2]
group1_samples <- metadata$Sample_ID[metadata$Condition == group1]
group2_samples <- metadata$Sample_ID[metadata$Condition == group2]
expression_subset <- gene_expression_matrix[, c(group1_samples, group2_samples)]

calculate_LFC <- function(x) {
  mean(x[group2_samples]) - mean(x[group1_samples])
}
log_fold_changes <- apply(expression_subset, 1, calculate_LFC)

factor_groups <- factor(c(rep(1, length(group1_samples)), rep(2, length(group2_samples))))
p_values <- rowttests(as.matrix(expression_subset), factor_groups)$p.value
adjusted_p_values <- mt.rawp2adjp(p_values, proc = "BH")$adjp[, 2]


results <- data.frame(
  Log_Fold_Change = log_fold_changes,
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values
)

# Ensure the column names match the metadata
colnames(metadata) <- c("Sample_ID", "Condition")

# Validate sample identifiers in metadata and gene expression matrix
valid_samples <- intersect(metadata$Sample_ID, colnames(gene_expression_matrix))
if (length(valid_samples) < length(metadata$Sample_ID)) {
  warning("Some samples in metadata are not present in the gene expression matrix. Only matching samples will be used.")
}

# Subset metadata and gene expression matrix
metadata <- metadata[metadata$Sample_ID %in% valid_samples, ]
gene_expression_matrix <- gene_expression_matrix[, valid_samples, drop = FALSE]

# Define groups
groups <- unique(metadata$Condition)
if (length(groups) < 2) stop("At least two conditions are required for differential expression analysis.")
group1 <- groups[1]
group2 <- groups[2]

# Identify samples for each group
group1_samples <- metadata$Sample_ID[metadata$Condition == group1]
group2_samples <- metadata$Sample_ID[metadata$Condition == group2]

# Subset the expression matrix
expression_subset <- gene_expression_matrix[, c(group1_samples, group2_samples), drop = FALSE]

# Calculate log-fold changes
calculate_LFC <- function(x) {
  mean(x[colnames(expression_subset) %in% group2_samples]) - mean(x[colnames(expression_subset) %in% group1_samples])
}
log_fold_changes <- apply(expression_subset, 1, calculate_LFC)

# Perform t-tests
factor_groups <- factor(c(rep(1, length(group1_samples)), rep(2, length(group2_samples))))
p_values <- rowttests(as.matrix(expression_subset), factor_groups)$p.value
adjusted_p_values <- mt.rawp2adjp(p_values, proc = "BH")$adjp[, 2]

# Compile results
results <- data.frame(
  Log_Fold_Change = log_fold_changes,
  P_Value = p_values,
  Adjusted_P_Value = adjusted_p_values
)

# Identify DEGs
significant_DEGs <- results[abs(results$Log_Fold_Change) > log2(2) & results$Adjusted_P_Value < 0.05, ]
write.table(rownames(significant_DEGs), file = "DEGs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Volcano plot
EnhancedVolcano(
  results,
  lab = rownames(results),
  x = 'Log_Fold_Change',
  y = 'Adjusted_P_Value',
  pCutoff = 0.05,
  FCcutoff = 2
)

# Heatmap for top 100 DEGs

#top_100_DEGs <- head(rownames(significant_DEGs), 100)
heatmap_data <- gene_expression_matrix[top_100_DEGs, ]
Heatmap(
  heatmap_data,
  name = "Expression",
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 10),
  top_annotation = HeatmapAnnotation(Condition = metadata$Condition)
)

# PCA for top 100 DEGs
top_100_pca <- prcomp(t(heatmap_data), center = TRUE, scale. = TRUE)
autoplot(top_100_pca, data = metadata, colour = 'Condition')

# Load Required Libraries
install.packages("readxl")   # For reading Excel files
install.packages("dplyr")    # For data manipulation

# Load the Libraries
library(readxl)

# Load datasets
sars_dataset <- read_excel("path_to_file/SARs dataset.xlsx")
alzheimers_dataset <- read_excel("path_to_file/Alzhaimer dataset.xlsx")

# Inspect the Data
head(sars_dataset)
head(alzheimers_dataset)

# Extract Gene Columns
# Replace 'GeneColumnName' with the actual column name
sars_genes <- sars_dataset$GeneColumnName
alzheimers_genes <- alzheimers_dataset$GeneColumnName

# Find Common Genes
common_genes <- intersect(sars_genes, alzheimers_genes)

# Save the Results
write.csv(common_genes, "common_genes.csv", row.names = FALSE)

# Verify the Results
common_genes_data <- read.csv("common_genes.csv")
head(common_genes_data)