
library(plyr)
library(data.table)
library(tidyverse)
library(fuzzyjoin)

###############################################################################
# Step 1. Specify File Locations
###############################################################################

### File locations ###

# I downloaded files locally, per sample, before plotting
main_dir = "C:/Users/addie/OneDrive - University of North Carolina at Chapel Hill/LAB/tardigrade_heterozygosity/aseq_results"
setwd(main_dir)

plots_dir = "C:/Users/addie/OneDrive - University of North Carolina at Chapel Hill/LAB/tardigrade_heterozygosity/aseq_results/PLOTS"


### ASEQ Results ###

## Identify variant files (one per RNA-seq sample)
ASEQ_variant_files <- dir(main_dir, recursive=T, include.dirs=T, pattern=".ASE.ASEQ", full.names=T)
ASEQ_variant_files

## Extract info from file path
samples <- c()
for (file in ASEQ_variant_files) {
  samples <- c(samples, unlist(strsplit(basename(file), ".ASE.ASEQ"))[1])
}
samples


### SnpEff-identified variants of interest ###

snpEff_variant_file <- "select_high_impact_variants.csv"


## Ranked ordered list of genes with high-impact variants based on expression level (from Plot_RNA_Expression_Data.R) ###

RNAseq_dir <- "C:/Users/addie/OneDrive - University of North Carolina at Chapel Hill/LAB/tardigrade_heterozygosity/rna_expression_analysis_results"
RNAseq_file <- paste0(RNAseq_dir, "/high_impact_variant_genes_ranked_by_expression.csv")


###############################################################################
# Step 2. Read in & Compile Data
###############################################################################

### PER-GENE ASEQ DATA ###

genes_stats <- read_tsv("genes_statistics.csv") # for all tested genes, fraction of SNPs with ASE, per sample
ase_genes   <- read_tsv("ASE_genes.csv") # list of genes with significant ASE, per sample
snps_per    <- read_tsv("snps_per_gene.csv") # number of SNPs per gene, per sample

genes_all <- genes_stats %>%
  left_join(ase_genes, by = c("chr", "start", "end", "gene")) %>%
  left_join(snps_per,  by = c("chr", "start", "end", "gene"),
            suffix = c("_fracASE", "_nSNP"))
head(genes_all)

genes_long <- genes_all %>%
  pivot_longer(
    cols = matches("SRR5218.*_(fracASE|nSNP)$"),
    names_to = c("sample", ".value"),
    names_pattern = "(SRR5218[0-9]+\\.sorted)_(fracASE|nSNP)"
  )
head(genes_long)

ase_counts_per_sample <- genes_long %>%
  mutate(ASE_call = fracASE > 0.5) %>%
  group_by(sample) %>%
  summarize(
    n_ASE_genes = sum(ASE_call),
    n_genes     = n(),
    prop_ASE    = n_ASE_genes / n_genes,
    .groups = "drop"
  )
ase_counts_per_sample

gene_intervals <- genes_all %>%
  distinct(chr, start, end, gene) %>%
  mutate(start = as.numeric(start),
         end   = as.numeric(end))
head(gene_intervals)


### PER-VARIANT ASEQ DATA ###

ASEQ_variants <- as.data.frame(NULL)
for (i in 1:length(ASEQ_variant_files)) {
  mydata <- read_tsv(ASEQ_variant_files[i])
  
  mydata$sample <- samples[i]
  mydata <- mydata %>%
    select(sample, chr, pos, ref, alt, A, C, G, T, af, RD, cov) %>%
    mutate(pos = as.numeric(pos))
  
  ASEQ_variants <- rbind(ASEQ_variants, mydata)
}
head(ASEQ_variants)

# Combine with gene intervals to identify which variants fall within which genes
ASEQ_variants <- bind_rows(
  lapply(unique(ASEQ_variants$chr), function(ch) {
    var_chr  <- filter(ASEQ_variants, chr == ch)
    gene_chr <- filter(gene_intervals, chr == ch)
    
    if (nrow(gene_chr) == 0) return(var_chr)  # no genes for this chr
    
    interval_inner_join(var_chr, gene_chr,
                        by = c("pos" = "start", "pos" = "end"),
                        maxgap = 0)
  })
) %>%
  select(sample, chr = chr.x, pos, gene, ref, alt, A, C, G, T, af, RD, cov)
head(ASEQ_variants)

# Calculate minor allele frequency based on read counts
ASEQ_variants <- ASEQ_variants %>%
  mutate(
    major_AF = pmax(A, C, G, T) / cov,
    minor_AF = 1 - major_AF
  )
head(ASEQ_variants)

# add column from ase_genes to indicate in how many samples a gene has significant ASE
ASEQ_variants <- ASEQ_variants %>%
  left_join(ase_genes %>% select(gene, samples.marked.as.ASE), by = "gene") %>%
  mutate(samples_marked_as_ASE = ifelse(is.na(samples.marked.as.ASE), 0, samples.marked.as.ASE)) %>%
  select(-samples.marked.as.ASE)
head(ASEQ_variants)


# ## generate version of ASEQ_variants where indels in REF or ALT allele are removed
# # remove rows where ALT or REF has greater than one letter
# ASEQ_variants_no_indels <- ASEQ_variants %>%
#   filter(nchar(ref) == 1 & nchar(alt) == 1)


### Read in SnpEff variants of interest ###

snpEff_variants = read.csv(snpEff_variant_file, header=T)
head(snpEff_variants)

# rename columns to match ASEQ files for merging
snpEff_variants <- snpEff_variants %>%
  rename(chr = scaffold, pos = position, ref = REF, alt = ALT, gene = Gene_ID) %>%
  select(chr, pos, ref, alt, gene, Annotation, Annotation_Impact)
head(snpEff_variants)


### Read in ranked list of genes with high-impact variants based on expression level ###

high_impact_variant_genes_ranked = read.csv(RNAseq_file, header=T)
head(high_impact_variant_genes_ranked)

#rename columns to match ASEQ files for merging
high_impact_variant_genes_ranked <- high_impact_variant_genes_ranked %>%
  rename(gene = gene_id) %>%
  select(gene, mean_TPM, genome_percentile_meanTPM)
head(high_impact_variant_genes_ranked)


### Combine ASEQ & RNA-seq info for final list of HIGH-IMPACT VARIANTS ###

ASEQ_variants_of_interest <- ASEQ_variants %>%
  # filter ASEQ list to only those variants that are also in the high-impact variants list
  inner_join(snpEff_variants, by = c("chr", "pos", "ref", "alt", "gene")) %>%
  # # add column to indicate whether each variant is in a gene with significant ASE (aka listed in ase_genes)
  # left_join(ase_genes %>% select(chr, start, end, gene) %>% distinct(), by = c("chr", "gene")) %>%
  # mutate(ASE_gene = ifelse(!is.na(start) & !is.na(end) & pos >= start & pos <= end, TRUE, FALSE)) %>%
  # select(-start, -end) %>%
  # add mean_TPM info for each gene
  left_join(high_impact_variant_genes_ranked, by = "gene")
head(ASEQ_variants_of_interest)


### Combine ASEQ & SnpEff & RNA-seq info for final list of VARIANTS WITHIN GENES CONTAINING HIGH-IMPACT VARIANTS ###
ASEQ_all_variants_within_genes_of_interest <- ASEQ_variants %>%
  # filter ASEQ list to only those genes that contain high-impact variants (aka listed in high_impact_variant_genes_ranked)
  inner_join(high_impact_variant_genes_ranked, by = "gene")
  # # add column to indicate whether each gene has significant ASE (aka listed in ase_genes)
  # left_join(ase_genes %>% select(chr, start, end, gene) %>% distinct(), by = c("chr", "gene")) %>%
  # mutate(
  #   ASE_gene = ifelse(!is.na(start) & !is.na(end) & pos >= start & pos <= end, TRUE, FALSE)) %>%
  # select(-start, -end)
head(ASEQ_all_variants_within_genes_of_interest)
  

### add column to ASEQ_variants to indicate whether gene is a gene of interest
ASEQ_variants <- ASEQ_variants %>%
  left_join(high_impact_variant_genes_ranked %>% mutate(Gene_of_Interest = TRUE) %>% select(gene, Gene_of_Interest), by = "gene") %>%
  mutate(Gene_of_Interest = ifelse(is.na(Gene_of_Interest), FALSE, Gene_of_Interest))
head(ASEQ_variants)

###############################################################################
### Adjust settings for ALL plots ###
###############################################################################

# theme settings
theme_set(theme_classic())
theme_update(text=element_text(size=10),
             strip.background = element_blank(), # remove box around facet labels
             strip.text = element_text(face = "bold"), # and bold text
             legend.key.size = unit(0.1, "in") # adjust legend sizing
)

###############################################################################
# Step 3. Plots
###############################################################################

### Are genes with high-impact variants more likely to show allele-specific expression (ASE) than all tested genes? ###

# compare genes in ASEQ_all_variants_within_genes_of_interest vs all genes tested by ASEQ (genes in ASEQ_variants)
# bar plot of percentage of genes in each list that have significant ASE




### How many genes in final show significant ASE, or not?

ggplot(ASEQ_variants_of_interest)+
  geom_bar(aes(x = ASE_gene, fill = sample), position = "dodge")+
  # label bars with counts
  geom_text(stat='count', aes(x = ASE_gene, group = sample, label=..count..), 
            position = position_dodge(width = 0.9), vjust=-0.5)+
  labs(x = "Gene with Significant ASE?", y = "Number of Variants")


### Are genes with high-impact variants more likely to show allele-specific expression (ASE)? ###

# compare ASEQ_variants vs ASEQ_all_variants_within_genes_of_interest
# bar plot of percentage of genes in each list that have significant ASE


### of the genes in this list with significant ASE, which allele is more highly expressed? (af = frequency of ALT allele)

## histogram of minor allele freq per variant/gene (mean across samples)
plotting_data <- ASEQ_variants_of_interest %>%
  filter(ASE_gene == TRUE) %>%
  # calculate mean af across samples for each chr/pos/gene
  group_by(chr, pos, gene) %>%
  summarize(mean_minor_AF = mean(minor_AF),
            sd_maf = sd(minor_AF),
            se_maf = sd_maf / sqrt(n()),
            .groups = "drop")
head(plotting_data)
length(unique(plotting_data$gene))
ggplot(plotting_data)+
  geom_histogram(aes(x = mean_minor_AF), bins = 40)+
  scale_x_continuous(limits = c(0,1))+
  labs(x = "Mean Minor Allele Frequency (across samples)", y = "Number of Variants")

## per variant and per sample, plot MAF of genes with statistically significant ASE
plotting_data <- ASEQ_variants_of_interest %>%
  filter(ASE_gene == TRUE)
# label samples as 1,2,3 instead of SRR5218xxx
plotting_data$sample <- factor(plotting_data$sample,
                               levels = unique(plotting_data$sample),
                               labels = c("1", "2", "3"))
head(plotting_data)
ggplot(plotting_data)+
  geom_point(aes(x = sample, y = minor_AF, color = sample), position = position_jitter(width = 0.2, height = 0))+
  facet_wrap(~gene)+
  scale_y_continuous(limits = c(0,0.5))+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
  labs(x = "Sample", y = "Minor Allele Frequency")


## per variant and sample, plot MAF of ALL alleles in final list
plotting_data <- ASEQ_variants_of_interest %>%
  # label samples as 1,2,3 instead of SRR5218xxx
  mutate(sample = factor(sample,
                         levels = unique(sample),
                         labels = c("1", "2", "3"))) %>%
  # re-order genes by mean af across samples
  group_by(chr, pos, gene) %>%
  mutate(mean_af = mean(af)) %>%
  ungroup() %>%
  mutate(gene = reorder(gene, minor_AF))
head(plotting_data)
ggplot(plotting_data)+
  geom_point(aes(x = sample, y = minor_AF, color = sample), position = position_jitter(width = 0.2, height = 0))+
  facet_wrap(~gene)+
  scale_y_continuous(limits = c(0,0.5))+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
  # add a * to plots with significant ASE
  geom_text(data = plotting_data %>% filter(ASE_gene == TRUE) %>% distinct(gene),
            aes(x = 2, y = 0.4, label = "*"), color = "red", size = 5) +
  scale_color_manual(values = c("1" = "#CC79A7", "2" = "#009E73", "3" = "#56B4E9"))+
  labs(x = "Sample", y = "Minor Allele Frequency")


## per variant, plot mean MAF across samples +/- SE
plotting_data <- ASEQ_variants_of_interest %>%
  # calculate mean MAF across samples for each chr/pos/gene
  group_by(chr, pos, gene) %>%
  summarize(mean_minor_AF = mean(minor_AF),
            sd_maf = sd(minor_AF),
            se_maf = sd_maf / sqrt(n()),
            .groups = "drop") %>%
  # add back in column to indicate whether each gene has significant ASE
  left_join(ase_genes %>% select(start, end, gene) %>% distinct(), by = "gene") %>%
  mutate(ASE_gene = ifelse(!is.na(start) & !is.na(end) & pos >= start & pos <= end, TRUE, FALSE)) %>%
  select(-start, -end) %>%
  # order genes by mean af
  mutate(gene = reorder(gene, mean_minor_AF))
head(plotting_data)
ggplot(plotting_data)+
  geom_point(aes(x = gene, y = mean_minor_AF))+
  geom_errorbar(aes(x = gene, ymin = mean_minor_AF - se_maf, ymax = mean_minor_AF + se_maf), width = 0.2)+
  scale_y_continuous(limits = c(0,0.5))+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
  # add a * to genes with significant ASE
  geom_text(data = plotting_data %>% filter(ASE_gene == TRUE),
            aes(x = gene, y = mean_minor_AF+se_maf+0.01, label = "*"), color = "red", size = 5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x = "Gene", y = "Mean Minor Allele Frequency")


### How does minor allele freq / ASE vary with transcript abundance (mean_TPM)? ###

# ## order af vs gene plot by mean_TPM instead of mean_af
# plotting_data <- ASEQ_variants_of_interest %>%
#   # calculate mean af across samples for each chr/pos/gene
#   group_by(chr, pos, gene) %>%
#   summarize(mean_af = mean(af),
#             sd_af = sd(af),
#             se_af = sd_af / sqrt(n()),
#             .groups = "drop") %>%
#   # add back in column to indicate whether each gene has significant ASE
#   left_join(ase_genes %>% select(start, end, gene) %>% distinct(), by = "gene") %>%
#   mutate(ASE_gene = ifelse(!is.na(start) & !is.na(end) & pos >= start & pos <= end, TRUE, FALSE)) %>%
#   select(-start, -end) %>%
#   # add back in mean_TPM info for each gene
#   left_join(high_impact_variant_genes_ranked, by = "gene") %>%
#   # order genes by mean TPM
#   mutate(gene = reorder(gene, mean_TPM))
# ggplot(plotting_data)+
#   geom_point(aes(x = gene, y = mean_af))+
#   geom_errorbar(aes(x = gene, ymin = mean_af - se_af, ymax = mean_af + se_af), width = 0.2)+
#   scale_y_continuous(limits = c(0,1))+
#   geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
#   # add a * to genes with significant ASE
#   geom_text(data = plotting_data %>% filter(ASE_gene == TRUE),
#             aes(x = gene, y = mean_af+se_af+0.01, label = "*"), color = "red", size = 5) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
#   labs(x = "Gene (ordered by mean TPM)", y = "Mean ALT Allele Frequency (af)")


## scatterplot of mean minor allele freq per gene of interest vs mean_TPM per gene of interest (if af < 1, MAF = af, if af > 1, MAF = 1 - af)
plotting_data <- ASEQ_all_variants_within_genes_of_interest %>%
  # calculate mean minor allele freq across variants, then samples for each gene
  group_by(gene, sample) %>%
  summarize(per_sample_mean_MAF = mean(minor_AF),
            per_sample_sd_MAF = sd(minor_AF),
            per_sample_se_MAF = per_sample_sd_MAF / sqrt(n()),
            .groups = "drop") %>%
  group_by(gene) %>%
  summarize(mean_MAF = mean(per_sample_mean_MAF),
            sd_MAF = sd(per_sample_mean_MAF),
            se_MAF = sd_MAF / sqrt(n()),
            .groups = "drop") %>%
  # add back in column to indicate whether each gene has significant ASE (TRUE is gene is listed in ase_genes)
  # do not use start/end, only gene name, and keep both TRUE and FALSE genes
  left_join(ase_genes %>% mutate(ASE_gene = TRUE) %>% select(gene, ASE_gene), by = "gene") %>%
  mutate(ASE_gene = ifelse(is.na(ASE_gene), FALSE, ASE_gene)) %>%
  # add back in mean_TPM info for each gene
  left_join(high_impact_variant_genes_ranked, by = "gene")
head(plotting_data)
ggplot(plotting_data)+
  geom_point(aes(x = mean_TPM, y = mean_MAF, color = ASE_gene))+
  # trend line
  geom_smooth(aes(x = mean_TPM, y = mean_MAF), method = "lm", se = FALSE, color = "blue", alpha = 0.3)+
  geom_errorbar(aes(x = mean_TPM, ymin = mean_MAF - se_MAF, ymax = mean_MAF + se_MAF, color = ASE_gene), width = 0.2)+
  scale_y_continuous(limits = c(0,0.5))+
  scale_x_continuous(trans = "log10")+
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
  labs(x = "Mean TPM (transcript abundance)", y = "Mean Minor Allele Frequency (MAF)", color = "Significant ASE in >= 1 sample?")



### Compare minor allele frequencies of variants in genes with significant ASE vs those without significant ASE ###
plotting_data <- ASEQ_variants %>%
  # calculate mean minor allele freq across variants, then samples for each gene
  group_by(gene, sample) %>%
  summarize(per_sample_mean_MAF = mean(minor_AF),
            per_sample_sd_MAF = sd(minor_AF),
            per_sample_se_MAF = per_sample_sd_MAF / sqrt(n()),
            .groups = "drop") %>%
  group_by(gene) %>%
  summarize(mean_MAF = mean(per_sample_mean_MAF),
            sd_MAF = sd(per_sample_mean_MAF),
            se_MAF = sd_MAF / sqrt(n()),
            .groups = "drop") %>%
  # add back in column to indicate whether each gene has significant ASE (TRUE is gene is listed in ase_genes)
  # do not use start/end, only gene name, and keep both TRUE and FALSE genes
  left_join(ase_genes %>% mutate(ASE_gene = TRUE) %>% select(gene, ASE_gene), by = "gene") %>%
  mutate(ASE_gene = ifelse(is.na(ASE_gene), FALSE, ASE_gene))
head(plotting_data)
# histogram of mean MAF for genes with significant ASE vs those without
ggplot(plotting_data)+
  geom_histogram(aes(x = mean_MAF, fill = ASE_gene, color = ASE_gene), position = "identity", alpha = 0.5, bins = 40)+
  scale_x_continuous(limits = c(0,0.5))+
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black"), labels = c("TRUE" = "Significant ASE", "FALSE" = "Not Significant ASE"))+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), labels = c("TRUE" = "Significant ASE", "FALSE" = "Not Significant ASE"))+
  labs(x = "Mean Minor Allele Frequency (MAF)", y = "Number of Genes",
       fill = "Significant ASE in >= 1 sample?", color = "Significant ASE in >= 1 sample?")


### Compare minor allele frequencies of genes with high-impact variants with all genes tested by ASEQ ###
plotting_data <- ASEQ_variants %>%
  # calculate mean minor allele freq across variants, then samples for each gene
  group_by(gene, sample) %>%
  summarize(per_sample_mean_MAF = mean(minor_AF),
            per_sample_sd_MAF = sd(minor_AF),
            per_sample_se_MAF = per_sample_sd_MAF / sqrt(n()),
            .groups = "drop") %>%
  group_by(gene) %>%
  summarize(mean_MAF = mean(per_sample_mean_MAF),
            sd_MAF = sd(per_sample_mean_MAF),
            se_MAF = sd_MAF / sqrt(n()),
            .groups = "drop") %>%
  # add TRUE/FALSE column for whether gene contains high-impact variant of interest
  left_join(high_impact_variant_genes_ranked %>% mutate(High_Impact_Variant_Gene = TRUE) %>% select(gene, High_Impact_Variant_Gene), by = "gene") %>%
  mutate(High_Impact_Variant_Gene = ifelse(is.na(High_Impact_Variant_Gene), FALSE, High_Impact_Variant_Gene))
head(plotting_data)
# histogram of mean MAF for genes with high-impact variants vs all genes tested by ASEQ
ggplot(plotting_data)+
  geom_histogram(aes(x = mean_MAF, fill = High_Impact_Variant_Gene, color = High_Impact_Variant_Gene),
                 position = "identity", alpha = 0.5, bins = 40)+
  facet_wrap(~High_Impact_Variant_Gene, nrow = 2, scales = "free_y",
             labeller = labeller(High_Impact_Variant_Gene = c("TRUE" = "Genes with high-impact variants", "FALSE" = "All genes tested by ASEQ")))+
  scale_x_continuous(limits = c(0,0.5))+
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "black"), labels = c("TRUE" = "Genes with high-impact variants", "FALSE" = "All genes tested by ASEQ"))+
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black"), labels = c("TRUE" = "Genes with high-impact variants", "FALSE" = "All genes tested by ASEQ"))+
  labs(x = "Mean Minor Allele Frequency (MAF)", y = "Number of Genes",
       fill = "Gene with high-impact variant?", color = "Gene with high-impact variant?")


### Histogram of allele frequencies for high-impact variants vs all variants tested by ASEQ
plotting_data <- ASEQ_variants %>%
  # calculate mean minor allele freq across samples for each variant
  group_by(chr, pos, gene) %>%
  summarize(mean_MAF = mean(minor_AF),
            sd_MAF = sd(minor_AF),
            se_MAF = sd_MAF / sqrt(n()),
            .groups = "drop") %>%
  # add TRUE/FALSE column for whether variant has high impact
  left_join(snpEff_variants %>% mutate(High_Impact_Variant = TRUE) %>% select(chr, pos, gene, High_Impact_Variant), by = c("chr", "pos", "gene")) %>%
  mutate(High_Impact_Variant = ifelse(is.na(High_Impact_Variant), FALSE, High_Impact_Variant))
head(plotting_data)
ggplot(plotting_data)+
  geom_histogram(aes(x = mean_MAF, fill = High_Impact_Variant, color = High_Impact_Variant),
                 position = "identity", alpha = 0.5, bins = 40)+
  facet_wrap(~High_Impact_Variant, nrow = 2, scales = "free_y",
             labeller = labeller(High_Impact_Variant = c("TRUE" = "High-impact variants", "FALSE" = "All variants tested by ASEQ")))+
  scale_x_continuous(limits = c(0,0.5))+
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "black"), labels = c("TRUE" = "High-impact variants", "FALSE" = "All variants tested by ASEQ"))+
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "black"), labels = c("TRUE" = "High-impact variants", "FALSE" = "All variants tested by ASEQ"))+
  labs(x = "Mean Minor Allele Frequency (MAF)", y = "Number of Variants",
       fill = "High-impact variant?", color = "High-impact variant?")


### Look at specific genes

# add to list of all variants a column indicating how many samples show ASE (ase_genes$samples.marked.as.ASE or 0 if not included)
all_variants <- ASEQ_variants %>%
  left_join(ase_genes %>% select(gene, samples.marked.as.ASE), by = c("gene")) %>%
  mutate(samples_marked_as_ASE = ifelse(is.na(samples.marked.as.ASE), 0, samples.marked.as.ASE)) %>%
  select(-samples.marked.as.ASE)
head(all_variants)

# randomly select 5 genes from all_variants with significant ASE in all 3 samples, and 5 with significant ASE in 0 samples
set.seed(123) # for reproducibility
genes_to_plot <- all_variants %>%
  group_by(gene) %>%
  summarize(samples_marked_as_ASE = max(samples_marked_as_ASE), .groups = "drop") %>%
  filter(samples_marked_as_ASE == 3 | samples_marked_as_ASE == 0) %>%
  group_by(samples_marked_as_ASE) %>%
  sample_n(5) %>%
  pull(gene)
genes_to_plot

# per gene in genes_to_plot, plot af for each variant (chr/pos) and sample, colored by whether gene has significant ASE in all 3 samples vs 0 samples
plotting_data <- all_variants %>%
  filter(gene %in% genes_to_plot) %>%
  # label samples as 1,2,3 instead of SRR5218xxx
  mutate(sample = factor(sample,
                         levels = unique(sample),
                         labels = c("1", "2", "3"))) %>%
  # order genes by mean af across samples
  group_by(chr, pos, gene) %>%
  mutate(mean_MF = mean(minor_AF)) %>%
  ungroup() %>%
  # flag genes with significant ASE
  mutate(gene = ifelse(samples_marked_as_ASE == 3,
                       paste0(gene, " *"), paste0(gene))) %>%
  # re-order by number of samples with significant ASE
  mutate(gene = reorder(gene, samples_marked_as_ASE))
head(plotting_data)
ggplot(plotting_data)+
  geom_point(aes(x = sample, y = minor_AF, color = cov), position = position_jitter(width = 0.2, height = 0))+
  facet_wrap(~gene, nrow = 2)+
  scale_y_continuous(limits = c(0,0.5))+
  scale_color_viridis_c(trans = "log10",
                        limits = c(1,NA)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
  labs(x = "Sample", y = "Minor Allele Frequency", color = "coverage")


### Repeat above (look at specific genes) but from short list of ASEQ_all_variants_within_genes_of_interest

head(ASEQ_all_variants_within_genes_of_interest)

# organize data
all_variants_of_interest <- ASEQ_all_variants_within_genes_of_interest %>%
  # add column indicating how many samples show ASE
  left_join(ase_genes %>% select(gene, samples.marked.as.ASE), by = c("gene")) %>%
  mutate(samples_marked_as_ASE = ifelse(is.na(samples.marked.as.ASE), 0, samples.marked.as.ASE)) %>%
  select(-samples.marked.as.ASE) %>%
  # calculate minor allele freq
  mutate(MAF = ifelse(af <= 0.5, af, 1 - af)) %>%
  # add mean_TPM info for each gene
  left_join(high_impact_variant_genes_ranked, by = "gene")
head(all_variants_of_interest)

# randomly select 5 genes with significant ASE in all 3 samples, and 5 with significant ASE in 0 samples
set.seed(456) # for reproducibility
genes_to_plot <- all_variants_of_interest %>%
  group_by(gene) %>%
  summarize(samples_marked_as_ASE = max(samples_marked_as_ASE), .groups = "drop") %>%
  filter(samples_marked_as_ASE == 3 | samples_marked_as_ASE == 0) %>%
  group_by(samples_marked_as_ASE) %>%
  sample_n(4) %>%
  pull(gene)
genes_to_plot

# per gene in genes_to_plot, plot af for each variant (chr/pos) and sample, colored by coverage (cov) per variant
plotting_data <- all_variants_of_interest %>%
  filter(gene %in% genes_to_plot) %>%
  # label samples as 1,2,3 instead of SRR5218xxx
  mutate(sample = factor(sample,
                         levels = unique(sample),
                         labels = c("1", "2", "3"))) %>%
  # order genes by mean af across samples
  group_by(chr, pos, gene) %>%
  mutate(mean_af = mean(af)) %>%
  ungroup() %>%
  # flag genes with significant ASE
  mutate(gene = ifelse(samples_marked_as_ASE == 3,
                       paste0(gene, " *"), paste0(gene))) %>%
  # re-order by number of samples with significant ASE
  mutate(gene = reorder(gene, samples_marked_as_ASE))
head(plotting_data)
ggplot(plotting_data)+
  geom_point(aes(x = sample, y = MAF, color = cov), position = position_jitter(width = 0.2))+
  facet_wrap(~gene, nrow = 2)+
  scale_y_continuous(limits = c(0,0.5))+
  scale_color_viridis_c(trans = "log10",
                        limits = c(1,2000)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")+
  labs(x = "Sample", y = "Minor Allele Frequency", color = "coverage")


### make the above plots (minor allele freq per gene and sample) for all genes containing high-impact variants ###
plotting_data <- ASEQ_all_variants_within_genes_of_interest %>%
  # label samples as 1,2,3 instead of SRR5218xxx
  mutate(sample = factor(sample,
                         levels = unique(sample),
                         labels = c("1", "2", "3"))) %>%
  # flag genes with significant ASE
  mutate(gene = ifelse(samples_marked_as_ASE == 3,
                       paste0(gene, " *"), paste0(gene))) %>%
  # re-order by number of samples with significant ASE
  ungroup() %>%
  mutate(gene = reorder(gene, samples_marked_as_ASE))



