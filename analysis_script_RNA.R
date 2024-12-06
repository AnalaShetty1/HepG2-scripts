library(Seurat)
library(gdata)
library(tidyverse)
library(ggplot2)
library(scales)
library(viridis)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(clusterProfiler)
library(ggpubr)
library(org.Hs.eg.db)
library(SeuratWrappers)
library(monocle3)
library(stringr)

setwd("/home/data/301177")

# Load data
x <- read.table("Counts/subread_counts.txt", header=T, row.names=1, sep="\t", stringsAsFactors = F)
# Keep gene lengths as first column
cnt <- x[,5:ncol(x)]

# Rename cell names with updated IDs
renaming <- read.table('RNA_sample_renaming.txt', header=T, sep='\t')
new_cols = c('Length')

for (c in seq(2,ncol(cnt))){
  new_col = colnames(cnt)[c]
  for (n in seq(1,nrow(renaming))){
    if (str_detect(colnames(cnt)[c], renaming$Old.Name[n])){
      new_col = str_replace(colnames(cnt)[c], renaming$Old.Name[n],renaming$New.Name[n])
      print(new_col)
    }
  }
  new_cols <- c(new_cols, new_col)
}

colnames(cnt) <- new_cols

# Gene length normalization and log1p transform
cnt.norm = cnt[,2:ncol(cnt)] / cnt[,1]
cnt.norm = t(t(cnt.norm) * 1e6 / colSums(cnt.norm))
cnt.norm.1p <- log1p(cnt.norm)

# Remove genes with 0 counts across samples
keep <- rowSums(cnt.norm) > 0

# Remove gene length column in count data and subset genes
rnacnt <- cnt[keep, 2:ncol(cnt)]
cnt.norm.1p <- cnt.norm.1p[keep, ]

genes <- read.table("Counts/gene_id_gene_name_map.txt", stringsAsFactors = FALSE)
names(genes) <- c("ensembl_gene_id", "mgi_symbol")
mitochondrial_genes <- genes[substr(genes$mgi_symbol, 1, 3) == 'MT-', ]$ensembl_gene_id

# if needed remove mitochondrial genes
#rnacnt = rnacnt[!rownames(rnacnt) %in% mitochondrial_genes,]
#cnt.norm.1p = cnt.norm.1p[!rownames(cnt.norm.1p) %in% mitochondrial_genes,]

# Create Seurat object with count data
sce <- CreateSeuratObject(counts=rnacnt, names.field=1)

# Add normalized slot
sce <- SetAssayData(object=sce, slot="data", new.data=cnt.norm.1p)

# Rename metadata
sce$nCount = sce$nCount_RNA
sce$nGene = sce$nFeature_RNA   
sce$mito_ratio <-  PercentageFeatureSet(object = sce, features = mitochondrial_genes) / 100

# Add phase information
phase <- sapply((strsplit(colnames(rnacnt), split = "_")), '[', 4)
phase <- ifelse(is.na(phase), 'S', phase)
sce$phase = phase  

# Filter cells with low genes / UMIs
sce <- subset(
  x = sce,
  subset =
    (nCount >= 1000) &
    (nGene >= 500)
)

# Process data with the standard Seurat workflow
sce <- FindVariableFeatures(object = sce)
sce <- ScaleData(object = sce)
sce <- RunPCA(object = sce, npcs = 43)

# Compute maximal number of PCs
pct <- sce[["pca"]]@stdev / sum(sce[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
max_pc <- min(co1, co2)

# UMAP, find neighbors and clustering with that number of PCs
sce <- RunUMAP(object = sce, dims = 1:max_pc)
sce <- FindNeighbors(object = sce)
sce <- FindClusters(object = sce)

# Plotting the PCA using the 100 most variable features
pca_plot <- DimPlot(sce,
                    reduction = 'pca', 
                    group.by = 'phase',
                    pt.size = 3) + theme_minimal() +
  scale_color_brewer(palette="Paired") +
  labs(title='PCA plot')

ggsave('output/pca_100.png', plot = pca_plot, dpi = 400 , bg = "white")
ggsave('output/pca_100.pdf', plot = pca_plot, dpi = 400 , bg = "white")

# Cell cycle markers
s_genes = c("MCM2","MCM6","ORC1","ORC2","CDC6","CCNE1","E2F1","E2F2","PCNA","RPA1","RPA2","RPA3")

# Get Ensembl IDs
ens_s = c()
for (gene in s_genes){
  if (gene %in% genes$mgi_symbol){
    ens = genes[genes$mgi_symbol==gene,]$ensembl_gene_id
    if (ens %in% rownames(sce)){
      ens_s <- c(ens_s, ens)
    }
  }
}

# G1/G2 markers are loaded with Seurat
ens_g = c()
for (gene in g2m.genes){
  if (gene %in% genes$mgi_symbol){
    ens = genes[genes$mgi_symbol==gene,]$ensembl_gene_id
    ens_g <- c(ens_g, ens)
  }
}

sce <- ScaleData(object = sce, features = rownames(sce))

# Computing sum count of markers
sce$sCount <- colSums(sce[c(ens_s),]@assays$RNA@scale.data)
sce$gCount <- colSums(sce[c(ens_g),]@assays$RNA@scale.data)

# S-phase markers plot
plot <- sce@meta.data %>%
  ggplot(aes(y=sCount, x=phase, fill=phase, group=phase)) + 
  geom_boxplot() + 
  theme_light() +
  xlab("Stage") +
  ylab("Scaled TPM normalized read counts") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  theme(axis.title.x =  element_text(vjust = -2)) +
  theme(axis.title.y =  element_text(vjust = 2)) +
  theme(plot.margin=unit(c(2,1,1,1), 'cm')) +
  ggtitle("S - phase marker genes expression") +
  scale_fill_brewer(palette='Paired') +
  stat_compare_means(label = "p.signif", method='t.test', comparisons=list(c('G1','G2'), c('G2','S'), c('G1','S')))

ggsave('output/s_gene_markers.png', plot = plot, dpi = 400 , bg = "white")
ggsave('output/s_gene_markers.pdf', plot)

# G1/G2 markers plot
plot <- sce@meta.data %>%
  ggplot(aes(y=gCount, x=phase, fill=phase, group=phase)) + 
  geom_boxplot() + 
  theme_light() +
  xlab("Stage") +
  ylab("Scaled TPM normalized read counts") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  theme(axis.title.x =  element_text(vjust = -2)) +
  theme(axis.title.y =  element_text(vjust = 2)) +
  theme(plot.margin=unit(c(2,1,1,1), 'cm')) +
  ggtitle("G1/G2 - phase marker genes expression") +
  scale_fill_brewer(palette='Paired') +
  stat_compare_means(label = "p.signif", method='t.test', comparisons=list(c('G1','G2'), c('G2','S'), c('G1','S')))

ggsave('output/g2m_gene_markers.png', plot = plot, dpi = 400 , bg = "white")
ggsave('output/g2m_gene_markers.png', plot)

# DE, GO enrichment, GSEA heatmaps between groups ##############################
# Set object idents
Idents(sce) = sce$phase

# Change the pairs for other groups -> G1 vs S and G2 vs S
cells.subset <- subset(sce, subset = (phase %in% c('G1', 'G2')))

# Find differentially expressed genes between groups
s_vs_g1 <- FindMarkers(cells.subset, ident.1 = 'G1', ident.2 = 'G2', test.use = 'DESeq2')
significant <- s_vs_g1[(!is.na(s_vs_g1$p_val_adj) & s_vs_g1$p_val_adj < 0.05),]

# Select the scaled data for the DE genes
s_vs_g1 = s_vs_g1[rownames(s_vs_g1) %in% rownames(cells.subset@assays$RNA@scale.data), ]

# Plot scaled top 100 DE genes between grops
heatmap <- DoHeatmap(cells.subset, features = rownames(s_vs_g1[1:100,]), draw.lines = F) +
  theme(plot.margin = margin(1,1,1,1, "cm")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  #scale_fill_viridis(na.value='white') +
  guides(colour=FALSE)

ggsave('output/G1-G2_heatmap.png', heatmap)
ggsave('output/G1-G2_heatmap.pdf', heatmap)

# GO-term enrichment on the significant terms
r_go <- clusterProfiler::enrichGO(
  gene = rownames(significant),
  universe = rownames(sce),
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = 'ALL',
  pAdjustMethod = "BH",
  readable = TRUE
)

write.table(r_go, 'output/GO_G1_vs_S.tsv')
ggsave('output/GO_G1_vs_S.png', barplot(r_go))
ggsave('output/GO_G1_vs_S.pdf', barplot(r_go))

# Load hallmark gene sets from msidbr
msigdbr_data <- msigdbr::msigdbr(species = 'human', category = 'H') %>%
  dplyr::mutate(gs_name = stringr::str_replace(gs_name, "HALLMARK_", ""))

# Reformat input data
colnames(s_vs_g1) <- c('p_val', 'log2FoldChange', 'pct.1', 'pct.2', 'padj')
deseq_data <- s_vs_g1 %>%
  dplyr::filter(!is.na(.data$log2FoldChange)) %>%
  dplyr::arrange(dplyr::desc(.data$log2FoldChange)) %>%
  dplyr::select(.data$log2FoldChange)

# Convert to vector
gsea_vector <- deseq_data$log2FoldChange %>% purrr::as_vector()
names(gsea_vector) <- rownames(deseq_data)

# GSEA analysis
r_gsea <- clusterProfiler::GSEA(
  geneList = gsea_vector,
  eps = 0,
  seed = T,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(msigdbr_data, .data$gs_name, .data$ensembl_gene)
)

write.table(r_gsea, 'output/GSEA_G1_vs_S.tsv')
ggsave('output/GSEA_G1_vs_S.png', ridgeplot(r_gsea))
ggsave('output/GSEA_G1_vs_S.pdf', ridgeplot(r_gsea))


# Compare to bulk RNA-seq
# Load data
rna = read.table('BulkRepliseqdata/GSE90322_norm_counts_FPKM_GRCh38.p13_NCBI.tsv', header = 1, row.names = NULL)

# Map gene IDs to Ensembl
mapping = AnnotationDbi::mapIds(org.Hs.eg.db,
                                keys = as.character(rna$GeneID),
                                column = "ENSEMBL",
                                keytype ="ENTREZID")
rna$GeneID = mapping

# Remove genes with no Ensembl ID
rna = rna[!is.na(rna$GeneID),]

# Merge scRNA and bulk HepG2 data
scrna = as.data.frame(sce@assays$RNA@scale.data)
scrna = rownames_to_column(scrna, 'GeneID')
df = inner_join(scrna, rna)
rownames(df) = make.names(df$GeneID, unique = TRUE)
df$GeneID = NULL

# Set grouping vector and plot
group_vector = c(rep('HepG2', length(Idents(sce))), rep('Bulk_HepG2', 2))
heatmap = Heatmap(df,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(Group = group_vector),
        show_row_dend = FALSE,
        show_row_names = F, show_column_names = F,
        name='Expression')

ggsave('output/bulk_vs_scRNA.png', heatmap)
ggsave('output/bulk_vs_scRNA.pdf', heatmap)
