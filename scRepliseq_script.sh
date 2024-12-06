# Load environment
conda activate screplseq

# Run analysis using scripts Step 3-7:

# Step 3 
for file in ../kronos_output/BAM/*.bam ; do n=$(basename "$file" .bam) ; bash scripts/Step3_load_mapped_reads.sh $file HepG2_output/bin_fragment/ $n test/hg38-blacklist.v2.formated.bed test/hg38.chrom.size ; done

# Step 4
bash scripts/Step4_MAD_score_QC.sh HepG2_output/bin_fragment/bins/ HepG2_output/mad/ 'HepG2'

# Step 5.1
for file in HepG2_output/bin_fragment/bins/*gDNA_G1* ; do bash scripts/Step5_1_Check_G1_cells.sh $file HepG2_output/G1_karyotypes/ ; done

# Step 5.2
bash scripts/Step5_2_Merge_G1_cells.sh HepG2_output/merged_G1/merged_G1.Rdata HepG2_output/bin_fragment/fragment/*gDNA_G1*

# Step 6
for file in HepG2_output/bin_fragment/fragment/* ; do bash scripts/Step6_log2_median_RT_scores.sh $file HepG2_output/RT/ HepG2_output/merged_G1/merged_G1.Rdata test/hg38-blacklist.v2.formated.bed  test/hg38.chrom.size ; done

# Step 7
for file in HepG2_output/bin_fragment/bins/* ; do bash scripts/Step7_Binarization.sh $file HepG2_output/binarized/ HepG2_output/merged_G1/merged_G1.Rdata test/hg38.chrom.size 200000 ; done
