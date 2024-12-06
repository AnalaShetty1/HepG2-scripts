library(Kronos.scRT)

# Mapping
Kronos.scRT::FastqToBam(
  bowtie2_index = '/home/data/bowtie-index/genome',
  File1 = list.files(path='trimmed/', pattern="_R1.fastq"),
  File2 = list.files(path='trimmed/', pattern="_R2.fastq"),
  outputdir = '/home/data/kronos_output',
  trim = F,
  cores = 6)

# Alignment statistics collection
bm_all = data.frame()
for (file in list.files(path='BAM/', pattern='.bam$')){
  print(file)
  bm <- Kronos.scRT::BamMetrics(paste0('BAM/', file), isPE = T)
  bm_all <- rbind(bm_all, bm)
}

# This table can be loaded in the next step
write.table(bm_all,file=paste0("alignment_stats/full_stats.tsv"), row.names = FALSE, quote=FALSE, sep='\t')

# Create bins from genome
bins_human = binning(
  RefGenome = '/home/data/bowtie-index/homo_sapiens_GRCh38_ensembl_release_107_ERCC_SIRV.fa',
  bowtie2_index = '/home/data/bowtie-index/genome',
  directory_to_bamfiles = '/home/data/kronos_output/BAM',
  cores = 10,
  bin_size = 20000
)

# Download and subset chromsizes
Chromsize = readr::read_tsv(url('http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes'), col_names = c('chr', 'size'))
Chromsize$chr = str_replace(Chromsize$chr, 'chr', '')
Chromsize = Chromsize[Chromsize$chr %in% seq(1,22),]

# CNV calling
SingleCell_human = Kronos.scRT::CallCNV(
  directory = '/home/data/kronos_output/BAM',
  chrom_size = Chromsize,
  bins = bins_human,
  basename = '301177',
  chr_prefix = '',
  chr_range=1:22,
  ploidy = 2,
  cores = 8
)

# Assign known cell phases
whoiswho = data.frame(Cell=SingleCell_human$PerCell$Cell, S_Phase = as.logical(str_detect(SingleCell_human$PerCell$Cell, 'gDNA_S')))
SingleCell_human$PerCell = Kronos.scRT::WhoIsWho(PerCell = SingleCell_human$PerCell, WhoIsWho = whoiswho)

# Open diagnostic window (previously assigned cell phases are not modified)
Diagnostic_output_human = Kronos.scRT::diagnostic(SingleCell_human$PerCell)

# Save diagnostic plots
ggsave('../plots/all_cells_plot.png',Diagnostic_output$all_cells_plot)
ggsave('../plots/first_filtering_plot.png',Diagnostic_output$first_filtering_plot)
ggsave('../plots/selected_G1_S_cells_plot.png',Diagnostic_output$selected_G1_S_cells_plot)
ggsave('../plots/S_phase_cell_distribution_plot.png',Diagnostic_output$S_phase_cell_distribution_plot)

# Print chosen settings for reference 
Diagnostic_output$Settings
#threshold_Sphase threshold_G1G2phase Sphase_first_part Sphase_second_part Ploidy RPMPH_TH RPM_TH basename group 
#<lgl>            <lgl>                           <dbl>              <dbl>  <dbl>    <int>  <dbl> <chr>    <chr> 
#  1 NA               NA                              0.907              0.713   2.27       75    170 301177   301177

# Adjust CNV values and filter cells based on chosen settings
SingleCell_human$PerCell = Kronos.scRT::AdjustPerCell(PerCell = SingleCell_human$PerCell,Settings = Diagnostic_output$Settings)
SingleCell_human$CNV = Kronos.scRT::AdjustCN(PerCell = SingleCell_human$PerCell, scCN = SingleCell_human$CNV)

# Genome bins generation
Bins_human = Kronos.scRT::GenomeBinning(
  Chr_size = Chromsize,
  Chr_filter =paste0('', 1:22),
  size = 20000,
  Cores = 6
)

# Rebin genome
SingleCell_human$SPhase = Kronos.scRT::Rebin(PerCell = SingleCell_human$PerCell, scCN = SingleCell_human$CNV, Bins = Bins_human, Sphase = T)
SingleCell_human$G1G2 = Kronos.scRT::Rebin(PerCell = SingleCell_human$PerCell, scCN = SingleCell_human$CNV, Bins = Bins_human, Sphase = F)
SingleCell_human$MedianG1G2=Kronos.scRT::BackGround(G1_scCN = SingleCell_human$G1G2)

# Replication state calculation
SingleCell_human$SPhase = Kronos.scRT::Replication_state(
  Samples = SingleCell_human$SPhase,
  background = SingleCell_human$MedianG1G2,
  Chr_filter = paste0('', 1:22),
  cores = 6
)

SingleCell_human$G1G2 = Kronos.scRT::Replication_state(
  Samples = SingleCell_human$G1G2,
  background = SingleCell_human$MedianG1G2,
  Chr_filter = paste0('', 1:22),
  cores = 6
)

# Pseudobulk computation
SingleCell_human$pseudobulk=Kronos.scRT::pseudoBulkRT(S_scCN = SingleCell_human$SPhase)

# Smooth RT
SingleCell_human$pseudobulk = SingleCell_human$pseudobulk %>% 
  group_by(chr) %>%
  mutate(RT=rollapply(RT, 4, mean, fill=0)) %>%
  ungroup()

# Generation of plots for all chromosomes
for (c in 1:nrow(Chromsize)){
  
  print(c)
  plot <- scRTplot(SingleCell_human$pseudobulk,
                   S_scCN = SingleCell_human$SPhase,
                   Coordinates = list(chr = Chromsize$chr[c], start = 0, end = Chromsize$size[c]),
                   rasterized_heatmap = T)
  
  ggsave(paste0('../plots/gDNA/chromosomes/',Chromsize$chr[c], '.png'), plot)
}

# Extraction of CNV per gene
CNV_annotated = Kronos.scRT::TW_GenomeAnnotation(SingleCell_human$CNV,
                                                 GenomeAnnotation = annotation) %>%
  filter(category != '_Unknown_') %>%
  select(category, copy_number, copy_number_corrected) %>% 
  group_by(category) %>%
  summarise(copy_number = mean(copy_number), copy_number_corrected = mean(copy_number_corrected)) %>%
  arrange(desc(copy_number_corrected))


