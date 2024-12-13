# Commands Used in Analysis
## Linux and Command Line Commands (Acquiring and Downloading of Data)

`#!/usr/bin/env bash`
`wget "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR2103774" -O SRR2103774.fastq.gz`

>*Note: commands were used for both datasets, unspecified here, majority of datasets were downloaded using `wget`, while those that were too large to download using that method were downloaded using the sra toolkit and the `fasterq-dump`*
`fastqc SRR2103774.fastq.gz`
>*Note: fastqc was run on all files individually*
`mamba install fastp`
`fastp -i SRR2103774.fastq.gz -o SRR2103774_trimmed.fastq.gz -h SRR2103774_trimmed_report.html -j SRR2103774_trimmed_report.json`
*Note: fastp summaries were checked before moving on*
`for fastq_file in *fastq.gz; do
sample=$(echo $fastq_file | sed 's/_fastq.gz//')
fastp -i "$fastq_file" \
-o "${sample}_trimmed.fastq.gz" \
-h "${sample}_trimmed_report.html" \
-j "${sample}_trimmed_report.json"
done`
`wget "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz" -O dm6.fa.gz
`
`wget "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/md5sum.txt" -O md5sum.txt`
`md5sum -c md5sum.txt`
`bwa index dm6.fa.gz`
`for file in *_trimmed.fastq.gz; do sample=$(echo $file | sed 's/_trimmed.fastq.gz//')
bwa mem dmelrel6_filtered.fa $file > ${sample}_trimmed.sam
done`
`mamba install subread`
`wget "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz" -O dmel-all-r6.48.gtf.gz`
`wget "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.48_FB2022_05/gtf/md5sum.txt" -O md5sum.txt`
`md5sum -c md5sum.txt`
`md5sum *.gz > md5check.txt`
`less md5check.txt`
`gunzip dmel-all-r6.48.gtf.gz`
`head dmel-all-r6.48.gtf`
`for file in *_trimmed.sam; do sample=$(echo $file | sed 's/_trimmed.samm//'); featureCounts -a dmel-all-r6.48.gtf -o ${sample}_counts.txt -t exon -g gene_id $file; done
`
`for file in *_counts.txt; do sample=$(echo $file | sed 's/_counts.txt//'); grep -v "^#" $file > "${sample}_cleaned_counts.txt"; done
`
`for file in *_cleaned_counts.txt; do
sample=$(echo $file | sed 's/_cleaned_counts.txt//')
bioawk -F'\t' '{n = split($7, counts, ";"); print $1, counts[1]}' "$file" > "${sample}_condensed_counts.txt"
done`

## R Analysis (Statistical Analysis and Plotting)
### Mattila Dataset
`setwd("~/Desktop/mattila_counts")`
`file_names <- list.files(pattern = "_condensed_counts.txt")`
`mattila_count_data_list <- list()`
`for (file in file_names) {temp_data <- read.delim(file, sep = " ", header = TRUE)  colnames(temp_data)[2] <- "Counts" dataset_name <- gsub("_condensed_counts.txt", "", file) mattila_count_data_list[[dataset_name]] <- temp_data }`
#### Preparing a count matrix
`rownames(mattila_count_matrix) <- mattila_count_matrix$Row.names`
`mattila_count_matrix <- mattila_count_matrix[, -1]`
`mattila_count_matrix <- NULL`
`for (dataset_name in names(mattila_count_data_list)) {
temp_data <- mattila_count_data_list[[dataset_name]]
if (!"Geneid" %in% colnames(temp_data)) {
stop("Geneid column not found in dataset: ", dataset_name)
}
rownames(temp_data) <- temp_data$Geneid
temp_data <- temp_data[, c("Geneid", "Counts"), drop = FALSE]
temp_Data <- temp_data[, "Counts", drop = FALSE]
colnames(temp_data) <- dataset_name
if (is.null(mattila_count_matrix)) {
mattila_count_matrix <- temp_data
} else {
mattila_count_matrix <- merge(mattila_count_matrix, temp_data, by = "row.names", all = TRUE)
rownames(mattila_count_matrix) <- mattila_count_matrix$Row.names
mattila_count_matrix <- mattila_count_matrix[, -1]
}
print(paste("Processed dataset:", dataset_name))
print(head(mattila_count_matrix))
}`
#### Preparing a metadata table
```{r}
mattila_metadata <- data.frame(
SampleName = c("LSD_WT_repA", "LSD_WT_repB", "LSD_WT_repC", "LSD_mlx1_repA", "LSD_mlx1_repB", "LSD_mlx1_repC", "HSD_WT_repA", "HSD_WT_repB", "HSD_WT_repC", "HSD_mlx1_repA", "HSD_mlxl1_repB", "HSD_mlx1_repC"),
Diet = c("LSD", "LSD", "LSD", "LSD", "LSD", "LSD", "HSD", "HSD", "HSD", "HSD", "HSD", "HSD"),
 Replicate = c("A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C"))
print(mattila_metadata)
```
#### Statistical Analysis
```{r}
BiocManager::install("DESeq2")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = mattila_count_matrix, colData = mattila_metadata, design = ~ Diet)
summary(res)
```
#### Changing gene identifiers to gene names
```{r}
url_gtf <- "https://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.48_FB2022_05/gtf/dmel-all-r6.48.gtf.gz"
destfile <- "dmel-all-r6.48.gtf.gz"
download.file(url_gtf, destfile, mode = "wb")
install.packages("R.utils")
library(R.utils)
gunzip("dmel-all-r6.48.gtf.gz")
gtf_data <- read.table("/Users/guanyinstevens/Desktop/mattila_counts/dmel-all-r6.48.gtf", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
gtf_data$gene_id <- gsub('.*gene_id ([^;]+).*', '\\1', gtf_data$V9)
gtf_data$gene_symbol <- gsub('.*gene_symbol ([^;]+).*', '\\1', gtf_data$V9)
head(gtf_data[, c("gene_id", "gene_symbol")])
gtf_cleaned <- gtf_data[, c("gene_id", "gene_symbol")]
lookup_table <- data.frame(gene_id = gtf_data$gene_id, gene_name = gtf_data$gene_symbol)
lookup_table <- data.frame(gene_id = gtf_data$gene_id, gene_name = gtf_data$gene_symbol)
View(lookup_table)
gene_id_column <- rownames(kokki_count_matrix)
gene_name_column <- lookup_table$gene_name[match(gene_id_column, lookup_table$gene_id)]
rownames(mattila_count_matrix) <- gene_name_column
head(mattila_count_matrix)
dds <- DESeqDataSetFromMatrix(countData = mattila_count_matrix, colData = mattila_metadata, design = ~ Diet)
dds <- DESeq(dds)
res <- results(dds)
summary(res)
```

#### Looking at HSD diets specifically. comparing by genotype
```{r}
hsd_metadata <- mattila_metadata[mattila_metadata$Diet == "HSD", ]
mlx1_samples <- hsd_metadata[grep("mlx1", hsd_metadata$SampleName, ignore.case = TRUE), ]
wt_samples <- hsd_metadata[grep("WT", hsd_metadata$SampleName, ignore.case = TRUE), ]
mlx1_samples <- hsd_metadata[grep("mlx1", hsd_metadata$SampleName, ignore.case = TRUE), ]
wt_samples <- hsd_metadata[grep("WT", hsd_metadata$SampleName, ignore.case = TRUE), ]
mlx1_counts <- mattila_count_matrix[, mlx1_samples$SampleName]
wt_counts <- mattila_count_matrix[, wt_samples$SampleName]
combined_HSD_counts <- cbind(mlx1_counts, wt_counts)
combined_HSD_metadata <- rbind(mlx1_samples, wt_samples)
combined_HSD_metadata$Genotype <- factor(
ifelse(grepl("mlx1", combined_HSD_metadata$SampleName), "mlx1", "WT"))
dds_HSD <- DESeqDataSetFromMatrix(countData = full_counts, colData = full_metadata, design = ~ Genotype)
full_counts <- cbind(mlx1_counts, wt_counts)
full_metadata <- rbind(mlx1_samples, wt_samples)
full_counts <- full_counts[, combined_HSD_metadata$SampleName]
dds_HSD <- DESeqDataSetFromMatrix(countData = full_counts, colData = combined_HSD_metadata, design = ~ Genotype)
dds_HSD <- DESeq(dds_HSD)
results(dds_HSD)
res_HSD <- results(dds_HSD)
significant_genes <- res_HSD[which(res$padj < 0.05), ]
View(significant_genes)
res_HSD.df <- as.data.frame(res_HSD)
res_HSD.df$gene_name <- rownames(res_HSD.df)
write_xlsx(res_HSD.df, "DESeq2_results_HSD_diet_Gene_Names.xlsx")
HSD_significant_genes.df <- as.data.frame(significant_genes)
write_xlsx(HSD_significant_genes.df, "Significant_Genes_HSD_Diet.xlsx")
HSD_significant_genes.df$gene_name <- rownames(HSD_significant_genes.df)
write_xlsx(HSD_significant_genes.df, "Significant_Genes_HSD_Diet.xlsx")
print(HSD_significant_genes.df)
```

#### Subsetting transcription factor genes from significant genes
```{r}
url <- "https://www.mrc-lmb.cam.ac.uk/genomes/FlyTF/all_candidates.csv"
destfile <- "FlyTF_all_candidates.csv"
download.file(url, destfile)
flytf_data <- read.table(destfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
symbols_only_df <- flytf_data[, "symbol", drop = FALSE]
print(symbols_only_df)
symbols_only_df$symbol <- tolower(trimws(symbols_only_df$symbol))
rownames(res) <- tolower(trimws(rownames(res)))
tf_names <- symbols_only_df$symbol
subset_HSD <- HSD_significant_genes.df[rownames(HSD_significant_genes.df) %in% tf_names, ]
```

### van den Berg/Kokki Dataset
> *Note: data acquisition, quality checks, trimming, alignment, and counts were done as explained above*
#### Creating count matrix and metadata table
```{r}
kokki_count_data_list <- list()
kokki_count_matrix <- NULL
for (file in file_list_1) {
temp_data <- read.delim(file, sep = "" "", header = TRUE)
colnames(temp_data)[2] <- "Counts"
dataset_name <- gsub("_condensed_counts.txt", " ", file)
kokki_count_data_list[[dataset_name]] <- temp_data
}
for (dataset_name in names(kokki_count_data_list)) {
temp_data <- kokki_count_data_list[[dataset_name]]
if (!"Geneid" %in% colnames(temp_data)) {
stop("Geneid column not found in dataset:", dataset_name)
}
temp_data <- temp_data[, c("Geneid", "Counts"), drop = FALSE]
rownames(temp_data) <- temp_data$Geneid
temp_data <- temp_data[, "Counts", drop = FALSE]  
colnames(temp_data) <- dataset_name
if (is.null(kokki_count_matrix)) {
kokki_count_matrix <- temp_data
} else {
kokki_count_matrix <- merge(kokki_count_matrix, temp_data, by = "row.names", all = TRUE)
rownames(kokki_count_matrix) <- kokki_count_matrix$Row.names
kokki_count_matrix <- kokki_count_matrix[, -1] 
}
print(paste("Processed dataset:", dataset_name))
}
kokki_metadata <- data.frame(
SampleName = c("kokki_WT_HSD_repA", "Kokki_WT_HSD_repB", "kokki_WT_HSD_repC", "Kokki_WT_HSD_repD", 
"kokki_WT_LSD_repA", "kokki_WT_LSD_repB", "Kokki_WT_LSD_repC", "Kokki_WT_LSD_repD"),
Diet = c("HSD", "HSD", "HSD", "HSD", "LSD", "LSD", "LSD", "LSD"),
Replicate = c("A", "B", "C", "D", "A", "B", "C", "D")
)
print(kokki_metadata)
``` 
#### Replacing gene identifiers with gene names
```{r}
gene_id_column <- rownames(kokki_count_matrix)
gene_name_column <- lookup_table$gene_name[match(gene_id_column, lookup_table$gene_id)]
rownames(kokki_count_matrix) <- gene_name_column
head(kokki_count_matrix)
```
#### Statistical analysis
```{r}
dds <- DESeqDataSetFromMatrix(countData = kokki_count_matrix, colData = kokki_metadata, design = ~ Diet)
res <- results(dds)
summary(res)
kokki_significant_genes <- res[which(res$padj < 0.05), ]
kokki_significant_genes.df$gene_name <- rownames(kokki_significant_genes)
kokki_significant_genes.df$gene_name <- rownames(kokki_significant_genes)
write_xlsx(kokki_significant_genes.df, "Kokki_Significant_Genes.xlsx")
kokki_top_genes <- head(rownames(log_kokki_transformed_significant_counts), 100)
kokki_top_genes_df <- as.data.frame((kokki_top_genes))
kokki_top_genes_df$gene_name <- rownames(kokki_top_genes)
write_xlsx(kokki_top_genes_df, "Kokki_top_genes.xlsx")
```

### Comparing the two datasets + identifying transcription factor genes
> *Note: transcription factors genes were identified from the same file as above*
```{r}
mattila_dataset <- read_excel("Significant_Genes_HSD_Diet.xlsx", sheet = "Sheet1")
kokki_dataset <- read_excel("Kokki_Significant_Genes.xlsx", sheet = "Sheet1")
merged_data <- merge(mattila_dataset, kokki_dataset, by = "gene_name")
write_xlsx(merged_data, "merged_data.xlsx")
significant_genes_both <- merged_data[!is.na(merged_data$padj.x) & merged_data$padj.x < 0.1 & !is.na(merged_data$padj.y) & merged_data$padj.y < 0.1, ]
significant_genes_both$gene_name <- tolower(significant_genes_both$gene_name)
tf_significant_genes_both <- significant_genes_both[significant_genes_both$gene_name %in% flytf_data$symbol, ]
heatmap_data <- tf_significant_genes_both[, c("log2FoldChange.x", "log2FoldChange.y")]
rownames(heatmap_data) <- tf_significant_genes_both$gene_name
```
#### Generating a heatmap for the significantly differentially expressed transcription factors
```{r}
pheatmap(heatmap_data,
scale = "none",
color = colorRampPalette(c("blue", "white", "red"))(100),
show_rownames = TRUE,
show_colnames = TRUE,
fontsize_row = 4,
main = "Heatmap of Log2 Fold Changes (No Scaling) for TFs Shared",
labels_col = c("Mattila Dataset", "Kokki Dataset")
)
```

