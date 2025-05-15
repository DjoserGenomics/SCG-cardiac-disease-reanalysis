# Loading R Packages
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(tibble)
library(biomaRt)
library(tximport)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(fgsea)
library(msigdbr)

# Creating Study Design file
studyDesign <- tibble(
  Accession = c("SRR24448340", "SRR24448339", "SRR24448338", "SRR24448337", "SRR24448336", "SRR24448335"),
  Sample = c("HS01", "HS02", "HS03", "CD01", "CD02", "CD03"),
  Condition = c("Healthy", "Healthy", "Healthy", "Disease", "Disease", "Disease")
)

# Ensuring Condition column in study design is categorical
studyDesign$Condition <- factor(studyDesign$Condition)

# Setting up Paths for kallisto abundance.tsv files
paths <- file.path('kallisto', studyDesign$Sample, 'abundance.tsv')
# Ensuring Paths are correct
all(file.exists(paths))

# Setting up Biomart to get Annotations 
myMart <- useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = 'https://oct2024.archive.ensembl.org')
mart.filters <- listFilters(myMart)

# Accessing v113 Human Annotations and selecting transcript id, gene id & gene name
annotations <- tryCatch({
  getBM(attributes = c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name'), mart = myMart) %>%
    as_tibble() %>%
    dplyr::rename(target_id = ensembl_transcript_id, gene_id = ensembl_gene_id, gene_name = external_gene_name)
}, error = function(e) {
  stop("Failed to fetch annotations from BioMart: ", e$message)
})

# Creating tx2gene that contains only transcript id and gene name for tximport
tx2gene <- annotations %>% dplyr::select(target_id, gene_name)

# Running tximport
txi <- tximport(paths,
                type='kallisto',
                tx2gene=tx2gene,
                txOut=FALSE,
                countsFromAbundance='lengthScaledTPM',
                ignoreTxVersion=TRUE)

# Removing Rows with '' or NA values
txi$counts <- txi$counts[!(rownames(txi$counts) == '' | is.na(rownames(txi$counts))), ]

# Keeping any gene where total counts across all samples is more than or equal to 10
keepers <- rowSums(txi$counts) >= 10
txi$counts <- txi$counts[keepers,]

# Assigning "Sample" column to be rownames for the study design
studyDesign <- studyDesign %>% column_to_rownames("Sample")

# Creating DESeq Dataset
dds <- DESeqDataSetFromTximport(txi = txi,
                              colData = studyDesign,
                              design = ~ Condition)
# Running DESeq
dds <- DESeq(dds)
# Producing Results and Setting Contrast for Comparison
res <- results(dds, contrast = c("Condition", "Disease", "Healthy"))

# Setting Cutoff variables
log2fcCutoff <- 1
padjCutoff <- 0.05

# Removing DEGs where their padj = NA
res.df <- res[!is.na(res$padj),] %>% as_tibble(rownames='gene_name')

# Filtering DEGs where abs(log2FoldChange) >= log2fcCutoff and padj < padjCutoff
DEGs <- res.df %>% 
  filter(padj < padjCutoff, abs(res.df$log2FoldChange) >= log2fcCutoff)

# Stabilizing Variance
vst <- vst(dds) 

# Functional Enrichment Analysis (FEA)

# GO Enrichment
ego <- enrichGO(gene = DEGs$gene_name, 
                OrgDb = "org.Hs.eg.db", 
                keyType = "SYMBOL", 
                ont = "ALL", 
                pvalueCutoff = 0.05, 
                qvalueCutoff = 0.05)

ego.df <- as_tibble(ego, rownames='OntologyID') %>% 
  arrange(p.adjust) %>% 
  filter(p.adjust < 0.05)

# Dotplot - top 8 terms
dotplot(ego, showCategory = 8, split = "ONTOLOGY", title = "Top GO Enrichment Terms: ~ BP/MF/CC")+ 
  facet_grid(ONTOLOGY ~ ., scales = "free")


gene.map <- bitr(DEGs$gene_name, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = "org.Hs.eg.db")

# KEGG Enrichment
ekg <- enrichKEGG(gene = gene.map$ENTREZID, 
                  organism = "hsa", 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05)

ekg.df <- as_tibble(ekg, rownames='hsa') %>% 
  arrange(p.adjust) %>% 
  filter(p.adjust < 0.05)

# Dotplot - top 10 pathways
dotplot(ekg, showCategory = 10, title = "Top 10 KEGG Pathways")


# GSEA
# Step 1: Prepare ranked gene list
ranks <- res.df$log2FoldChange
names(ranks) <- res.df$gene_name
ranks <- ranks[!is.na(ranks)] %>% sort(decreasing = TRUE)

# Step 2: Fetch MSigDB C2 KEGG gene sets
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "KEGG_LEGACY") %>%
  dplyr::select(gs_name, gene_symbol)
pathways <- split(hs_gsea_c2$gene_symbol, hs_gsea_c2$gs_name)

# Step 3: Run GSEA with fgsea
set.seed(123)
fgsea.res <- fgsea(pathways = pathways, 
                   stats = ranks, 
                   minSize = 15, 
                   maxSize = 500, 
                   nperm = 10000)

# Step 4: Process results
fgsea.df <- as_tibble(fgsea.res) %>%
  arrange(padj) %>%
  filter(padj < 0.05, abs(NES) > 1)

# KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION Enrichment Plot
top_pathway <- fgsea.df$pathway[2]
plotEnrichment(pathway = pathways[[top_pathway]], 
               stats = ranks) +
  labs(title = top_pathway, 
       subtitle = paste("NES =", round(fgsea.df$NES[2], 2), 
                        "p.adjust =", format(fgsea.df$padj[2], scientific = TRUE, digits = 3))) +
  theme_minimal()

# KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION Enrichment Plot
top_pathway <- fgsea.df$pathway[20]
plotEnrichment(pathway = pathways[[top_pathway]], 
               stats = ranks) +
  labs(title = top_pathway, 
       subtitle = paste("NES =", round(fgsea.df$NES[20], 2), 
                        "p.adjust =", format(fgsea.df$padj[20], scientific = TRUE, digits = 3))) +
  theme_minimal()

# Produce Plots
# Producing MA Plot (log2fc vs. mean)
plotMA(res, alpha = 0.05, main = "MA Plot: Disease vs Healthy (padj < 0.05)")

# Gathering PC Data & Producing a PCA Plot
pcaData <- plotPCA(vst, ntop=1000, "Condition", returnData = TRUE)

ggplot(pcaData, aes(PC1, PC2, color=Condition, label=name)) +
  geom_point(size=4) +  
  geom_text(aes(label = name), vjust = -1, size = 4, show.legend = FALSE) + 
  labs(x=paste0('PC1: ', round(100 * attr(pcaData, "percentVar")[1], 1), '%'),
       y=paste0('PC2: ', round(100 * attr(pcaData, "percentVar")[2], 1), '%'),
       title='PCA Plot with Sample Labels',
       subtitle='PC1 vs PC2') 

# Surrogate Variable Analysis (SVA) was tested to check for hidden batch effects
# However, PCA and clustering patterns remained unchanged after SVA correction,
# suggesting the structure is not due to technical variation.
# Code kept here for reproducibility, but not used in final DE analysis.

# mod <- model.matrix(~ condition, data=colData(dds))
# mod0 <- model.matrix(~ 1, data=colData(dds))
# svseq <- svaseq(as.matrix(counts(dds)), mod, mod0)
# colData(dds)$SV1 <- svseq$sv[, 1]
# design(dds) <- ~ SV1 + condition
# dds <- DESeq(dds)

# Producing a Volcano Plot
EnhancedVolcano(res.df,
                lab = res.df$gene_name, 
                x = "log2FoldChange", 
                y = "padj",
                pCutoff = 0.05, 
                FCcutoff = 1, 
                title = "Volcano Plot: Disease vs Healthy",
                pointSize = 2.0,
                xlim = c(-10, 10), 
                ylim = c(0, 10))

# Producing a Box plot showing per sample expression distribution
boxplot(log2(counts(dds, normalized=TRUE) + 1), main="Sample Expression Distribution")

# Producing Heatmaps
heatmapOrder <- c("HS01", "HS02", "HS03", "CD01", "CD02", "CD03")

# Inflammation Themed Heatmap

gseaInflammatoryData <- fgsea.df %>% 
  filter(fgsea.df$pathway == 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION') %>% 
  dplyr::select(leadingEdge)

gseaInflammatoryData <- gseaInflammatoryData$leadingEdge[[1]]

gseaInflammatoryGenes <- intersect(DEGs$gene_name, gseaInflammatoryData)

goInflammatoryData <- ego.df %>%
  filter(ego.df$ID == 'GO:0150076') %>%
  dplyr::select(geneID)

goInflammatoryGenes <- unlist(strsplit(goInflammatoryData$geneID[[1]], split="/"))

inflammatoryGenes <- unlist(unique(c(gseaInflammatoryGenes, goInflammatoryGenes)))

heatmapData <- assay(vst)[inflammatoryGenes, heatmapOrder]

pheatmap(heatmapData, cluster_cols = FALSE, main=paste0('Expression Patterns of Inflammation-Associated Genes'), show_rownames=TRUE)

# Adaptive Immunity Themed Heatmap
gseaImmunityData <- fgsea.df %>% 
  filter(fgsea.df$pathway == 'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION') %>% 
  dplyr::select(leadingEdge)

gseaImmunityData <- gseaImmunityData$leadingEdge[[1]]

gseaImmunityGenes <- intersect(DEGs$gene_name, gseaImmunityData)

goImmunityData <- ego.df %>%
  filter(ego.df$ID == 'GO:0019882') %>%
  dplyr::select(geneID)

goImmunityGenes <- unlist(strsplit(goImmunityData$geneID[[1]], split="/"))

immunityGenes <- unlist(unique(c(gseaImmunityGenes, goImmunityGenes)))

heatmapData <- assay(vst)[immunityGenes, heatmapOrder]

pheatmap(heatmapData, cluster_cols = FALSE, main=paste0('Expression Patterns of Adaptive Immunity Associated Genes'), show_rownames=TRUE)

# Producing a Sample-to-Sample distance heatmap
# Calculating Sample Distances
sampleDists <- dist(t(assay(vst)))

# Convert to matrix for visualization
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- colnames(vst)
rownames(sampleDistMatrix) <- colnames(vst)

pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         main = "Sample-to-Sample Distance Heatmap")
