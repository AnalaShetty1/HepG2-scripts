################### Part 2: scRepliseq analysis #####################################

# Load annotation
annotation = read.table('bowtie-index/annotation.bed')
colnames(annotation) <- c('chr', 'start', 'end', 'annotation')

# Initialize data frame
cell_scores = data.frame()

# Loop through scRepliseq output files
for (file in list.files('../scRepliseq-Pipeline/HepG2_output/RT', full.names = T)){
  cell = substr(basename(file),1,10)
  print(cell)
  
  # Load RT file and annotate
  rt = read.table(file)
  colnames(rt) = c('chr', 'start','end', 'RT')
  rt$chr = as.character(rt$chr)
  rt = Kronos.scRT::TW_GenomeAnnotation(Variability = rt, GenomeAnnotation = annotation)
  
  # Group by gene ID and rescale
  rt = rt %>% group_by(category) %>%
    filter(category != '_Unknown_') %>%
    summarise(RT = mean(RT)) %>% 
    mutate(RT=rescale(RT, to=c(0,1)))
  
  # Look for the corresponding sample in the RNA data
  rna_name = ''
  for (rna in colnames(cnt.norm)){
    if(grepl(cell, rna)){
      rna_name = rna
    }
  }
  
  if(rna_name != ''){
    
    # Extract RNA vector
    rna = data.frame('Expression' = cnt.norm[,rna_name])
    rna = rownames_to_column(rna, 'category')
    
    # Join RNA and RT vectors on gene ID
    df = inner_join(rna, as.data.frame(rt))
    
    # select genes with log2 FPKM > 1
    df$is_exp = log2(df$Expression) > 1
    
    # Bin genes in groups of 50 and compute
    # probability of expression
    df = df %>% arrange(RT) %>%
      group_by(x = ceiling(row_number()/50)) %>%
      summarise(prob = sum(is_exp)/50, RT=mean(RT))
    
    # Compute correlation line
    mod <- lm(df$prob ~ df$RT)
    cf <- coef(mod)
    
    # Save slope
    cell_scores = rbind(cell_scores, c(cell, cf[[2]], mean(df$RT)))
  }
}

# Save cell scores vector
colnames(cell_scores) = c('Cell', 'Slope', 'RT')
cell_scores$Slope = as.numeric(cell_scores$Slope)
cell_scores$RT = as.numeric(cell_scores$RT)
write.table(cell_scores, 'output/cell_scores.tsv', row.names = F, quote = F)

# Load pseudobulk RT
rt = read.table('../scRepliseq-Pipeline/HepG2_output/pseudobulk/RT/HepG2_pseudobulk_w200ks40k_map_count_median_log2.bedGraph')
colnames(rt) = c('chr', 'start','end', 'RT')
rt$chr = as.character(rt$chr)

# Annotate pseudobulk and group by gene ID
rt = Kronos.scRT::TW_GenomeAnnotation(Variability = rt, GenomeAnnotation = annotation)
rt = rt %>% group_by(category) %>%
  filter(category != '_Unknown_') %>%
  summarise(RT = mean(RT)) %>%
  mutate(RT=rescale(RT, to=c(0,1)))

# Binarize expression: FPKM > 1 in at least 1 sample
genes_expressed = data.frame(is_exp = rowSums(log2(cnt.norm) > 1) > 1,
                             exp = rowMeans(cnt.norm))
genes_expressed = rownames_to_column(genes_expressed, 'category')

# Join RNA and RT on gene ID
zyg = inner_join(genes_expressed, rt)

# Bin genes in groups of 50 and compute probability of expression
barplot_data <- zyg %>% arrange(RT) %>%
  group_by(x = ceiling(row_number()/50)) %>%
  summarise(prob = sum(is_exp)/50, RT=mean(RT))

# Compute correlation line
mod <- lm(barplot_data$prob ~ barplot_data$RT)
cf_pseudobulk <- coef(mod)

# Plot RT and probability of expression
rt_exp_plot <- ggplot(barplot_data, aes(x=RT, y=prob)) +
  geom_area(fill='darkgray') +
  geom_vline(xintercept = 0.5, linetype = 'dotted', linewidth=2, color='black', alpha=0.6) +
  geom_text(aes(x=0.7, y=0.8, label=paste0('Slope: ', as.character(round(cf_pseudobulk[[2]], 3))))) +
  ylab('Probability of expression per bin (log FPKM > 1 in at least one cell)') +
  xlab('Mean RT per bin (rounded)') +
  labs(title='HepG2') +
  theme_classic() +
  geom_smooth(method='lm', color='black', linewidth=0.8)

ggsave('output/expression_vs_RT_at_least_1_pseudobulk.png', rt_exp_plot)
ggsave('output/expression_vs_RT_at_least_1_pseudobulk.pdf', rt_exp_plot)

# Binarize RT into early and late
zyg$category <- ifelse(zyg$RT > 0.5, 'Early', 'Late')

# Generate boxplot for expression in early and late genes
rt_exp_boxplot <- ggplot(zyg, aes(x=category, y=exp, fill=category)) +
  geom_boxplot(outlier.size=2, weight=2) +
  scale_y_continuous(trans='log10') +
  scale_fill_brewer(palette = 'Set1') +
  ylab('Gene Expression (FPKM)') +
  xlab('RT category (early is > 0.5)') +
  labs(title='HepG2') +
  stat_compare_means(label.x.npc = 'center') 

ggsave('output/late_early_boxplot_at_least_1_pseudobulk.png', rt_exp_boxplot)
ggsave('output/late_early_boxplot_at_least_1_pseudobulk.pdf', rt_exp_boxplot)

# Load Bulk RT 
rt = read.table('BulkRepliseqdata/HepG2_R1_hg38_10kb_qNormSmoothed.bedgraph')
colnames(rt) = c('chr', 'start','end', 'RT')
rt$chr = gsub('chr', '', rt$chr)
rt$chr = as.character(rt$chr)

# Annotate RT and group by gene ID
rt = Kronos.scRT::TW_GenomeAnnotation(Variability = rt, GenomeAnnotation = annotation)
rt = rt %>% group_by(category) %>%
  filter(category != '_Unknown_') %>%
  summarise(RT = mean(RT)) %>%
  mutate(RT=rescale(RT, to=c(0,1)))

# Read in bulk RNA counts and map IDs to Ensembl
rna = read.table('BulkRepliseqdata/GSE90322_norm_counts_FPKM_GRCh38.p13_NCBI.tsv', header = 1)
mapping = AnnotationDbi::mapIds(org.Hs.eg.db,
                                keys = as.character(rna$GeneID),
                                column = "ENSEMBL",
                                keytype ="ENTREZID")

# Only keep genes that have Ensembl ID
rna$category = mapping
rna = rna[!is.na(rna$category),]

# Compute row means from the two samples
genes_expressed = data.frame(is_exp = rowSums(log2(rna[,2:3]) > 1) > 1,
                             exp = rowMeans(rna[,2:3]),
                             category = rna$category)

df = inner_join(genes_expressed, as.data.frame(rt))
barplot_data = df %>% arrange(RT) %>%
  group_by(x = ceiling(row_number()/50)) %>%
  summarise(prob = sum(is_exp)/50, RT=mean(RT))

mod_bulk <- lm(barplot_data$prob ~ barplot_data$RT)
cf_bulk <- coef(mod_bulk)

# Plot RT and probability of expression
rt_exp_plot <- ggplot(barplot_data, aes(x=RT, y=prob)) +
  geom_area(fill='darkgray') +
  geom_vline(xintercept = 0.5, linetype = 'dotted', linewidth=2, color='black', alpha=0.6) +
  geom_text(aes(x=0.7, y=0.8, label=paste0('Slope: ', as.character(round(cf_bulk[[2]], 3))))) +
  ylab('Probability of expression per bin (log FPKM > 1 in at least one cell)') +
  xlab('Mean RT per bin (rounded)') +
  labs(title='HepG2') +
  theme_classic() +
  geom_smooth(method='lm', color='black', linewidth=0.8)

ggsave('output/expression_vs_RT_at_least_1_bulk.png', rt_exp_plot)
ggsave('output/expression_vs_RT_at_least_1_bulk.pdf', rt_exp_plot)

# Binarize RT into early and late
df$category <- ifelse(df$RT > 0.5, 'Early', 'Late')

# Generate boxplot for expression in early and late genes
rt_exp_boxplot <- ggplot(df, aes(x=category, y=exp, fill=category)) +
  geom_boxplot(outlier.size=2, weight=2) +
  scale_y_continuous(trans='log10') +
  scale_fill_brewer(palette = 'Set1') +
  ylab('Gene Expression (FPKM)') +
  xlab('RT category (early is > 0.5)') +
  labs(title='HepG2') +
  stat_compare_means(label.x.npc = 'center') 

ggsave('output/late_early_boxplot_bulk.png', rt_exp_boxplot)
ggsave('output/late_early_boxplot_bulk.pdf', rt_exp_boxplot)

# Cell score plot
plot = ggplot(cell_scores, aes(x = 'HepG2', y=Slope, color=RT)) +
  geom_jitter(size=2.5) +
  geom_point(aes(y=cf_pseudobulk[[2]], shape="pseudobulk"), color='black', size=4, show.legend = T)+
  geom_point(aes(y=cf_bulk[[2]], shape="bulk"), color='black', size=3, show.legend = T)+
  xlab('') +
  scale_colour_gradient(low = 'red',high = 'blue') +
  labs(title = 'RT/RNA correlation line slope per cell') +
  scale_shape_manual(values = c(15, 17)) +
  labs(shape = "")

ggsave('output/correlation.png', plot)
ggsave('output/correlation.pdf', plot)

# Plot without bulk values
plot = ggplot(cell_scores, aes(x = 'HepG2', y=Slope, color=RT)) +
  geom_jitter(size=2.5) +
  xlab('') +
  scale_colour_gradient(low = 'red',high = 'blue') +
  labs(title = 'RT/RNA correlation line slope per cell')

ggsave('output/correlation_without_bulk.png', plot)
ggsave('output/correlation_without_bulk.pdf', plot)

# Plot slope vs RT
plot = ggplot(cell_scores, aes(x = RT, y=Slope, color=RT)) +
  geom_jitter(size=2.5) +
  xlab('RT') +
  scale_colour_gradient(low = 'red',high = 'blue') +
  labs(title = 'Correlation slope per cell vs RT') +
  geom_smooth(method='glm', color='black', linewidth=0.8)

ggsave('output/correlation_slope_vs_RT.png', plot)
ggsave('output/correlation_slope_vs_RT.pdf', plot)
