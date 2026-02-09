
library(plyr)
library(nlme)
library(scales)
library(ggplot2)
library(grid)
library(gridExtra)
library(data.table)
library(dplyr)
library(ggtext)

###############################################################################
# Specify File locations
###############################################################################

# I download the summary file locally before plotting
main_dir = "C:/Users/addie/OneDrive - University of North Carolina at Chapel Hill/LAB/tardigrade_heterozygosity/NCBI_genome_results"
setwd(main_dir)

plots_dir = "C:/Users/addie/OneDrive - University of North Carolina at Chapel Hill/LAB/tardigrade_heterozygosity/He_FIGURE_PLOTS"

###############################################################################
# Adjust settings for ALL plots
###############################################################################

# theme settings
theme_set(theme_classic())
theme_update(text=element_text(size=10),
             strip.background = element_blank(), # remove box around facet labels
             strip.text = element_text(face = "bold"), # and bold text
             legend.key.size = unit(0.1, "in") # adjust legend sizing
)

cbPalette <- c("#999999", "#56B4E9", "#E69F00", "#CC79A7", "#009E73")

###############################################################################
# Read in & Organize Data File(s)
###############################################################################

fileName <- "20251118_annotated_high_impact_variants_for_mRNA_search.csv"
# fileName <- "20251204_select_high_impact_variants_with_goslim.csv"
high_impact_variants = read.csv(fileName, header=T)
head(high_impact_variants)

# fileName <- "../Hypsibius_exemplaris_PredictedProteins.csv"
fileName <- "../Hypsibius_exemplaris_annotated_genes_with_goslim.csv"
ProteinAnnotations <- read.csv(fileName, header=T)
ProteinAnnotations$top_BLASTP_hit <- ifelse(ProteinAnnotations$top_BLASTP_hit == 0, NA, ProteinAnnotations$top_BLASTP_hit)
ProteinAnnotations$Uniprot_Protein_Annotation <- ifelse(ProteinAnnotations$Uniprot_Protein_Annotation=="Uncharacterized protein", NA, ProteinAnnotations$Uniprot_Protein_Annotation)
head(ProteinAnnotations)

fileName <- "../H_exemplaris_DNA_Repair_Genes_from_Clark-Hachtel_2024_tables_S3_S4.csv"
DNA_Repair_Genes <- read.csv(fileName, header=T)
head(DNA_Repair_Genes)

###############################################################################
# Analysis
###############################################################################

### Add protein annotations to my data ###
mydata <- high_impact_variants %>%
  left_join(ProteinAnnotations, by = "Gene_ID") %>%
  # left_join(DNA_Repair_Genes %>% select(Gene_ID, Protein, DNA_Repair_Pathway), by = "Gene_ID") %>%
  select(
    Gene_ID,
    Annotation,
    AA_pos,
    Ensembl_Description,
    top_BLASTP_hit,
    Uniprot_Protein_Annotation,
    Uniprot_subcellular_location,
    Gene_Ontology,
    GO_slim
    # Protein,
    # DNA_Repair_Pathway
  )
head(mydata)


## add column with simplified annotation classes
mydata <- mydata %>%
  mutate(
    Annotation_simple = case_when(
      grepl("stop_gained", Annotation) ~ "Stop gained",
      grepl("frameshift_variant", Annotation) ~ "Frameshift",
      grepl("splice_acceptor|splice_donor", Annotation) ~ "Splice site",
      grepl("start_lost|stop_lost", Annotation) ~ "Start/stop lost",
      TRUE ~ "Other high-impact"
    )
  )

### List of genes with predicted or known functions in any of the columns ###

descriptive_genes <- mydata %>%
  filter(!is.na(Ensembl_Description) | !is.na(top_BLASTP_hit) | !is.na(Uniprot_Protein_Annotation)) %>%
  select(Annotation, Gene_ID, Ensembl_Description, top_BLASTP_hit, Uniprot_Protein_Annotation)
head(descriptive_genes)
write.csv(descriptive_genes, "high_impact_variants_in_predicted_or_known_genes.csv", row.names = FALSE)
# Uniprot annotations fulfill the most genes, but sometimes a gene has a BLASTp hit but no Uniprot annotation


###############################################################################
### How many genes with >1 variant vs genes with only 1 variant? ###
###############################################################################

gene_counts <- mydata %>%
  count(Gene_ID, name = "n_rows")
head(gene_counts)

summary_counts <- gene_counts %>%
  count(n_rows, name = "n_genes")
summary_counts

ggplot(summary_counts)+
  geom_col(aes(x = as.factor(n_rows), y = n_genes))+
  geom_text(aes(label = n_genes, x = as.factor(n_rows), y = n_genes),
            vjust = -0.5) +
  # expand plot limits to fit text
  expand_limits(y = max(summary_counts$n_genes) * 1.1)+
  labs(x = "number of variants per gene", y = "number of genes")


### How many total variants? ###
n_variants <- nrow(mydata)
n_variants
### How many total genes? ###
n_genes <- n_distinct(mydata$Gene_ID)
n_genes

###############################################################################
### Gene counts with variant(s) per annotation category ###
###############################################################################

# prep data for plotting
gene_counts <- mydata %>%
  distinct(Gene_ID, Annotation_simple) %>%
  count(Annotation_simple, name = "n_genes")
variant_counts <- mydata %>%
  count(Annotation_simple, name = "n_variants")
plot_df <- left_join(gene_counts, variant_counts, by = "Annotation_simple")
# plot
ggplot(plot_df, aes(x = reorder(Annotation_simple, n_genes))) +
  # geom_col(aes(y = n_variants), fill = "grey85") +
  geom_col(aes(y = n_genes), fill = "grey40") +
  geom_text(
    aes(label = paste0(n_genes, " (", n_variants, ")"),
        y = n_genes),
    hjust = -0.1,
    size = 3
  ) +
  coord_flip() +
  expand_limits(y = max(plot_df$n_variants) * 1.2) +
  labs(
    x = NULL,
    y = "Numbers of genes"
    # title = "Genes with predicted high-impact heterozygous variants"
  )
plot_title <- "high_impact_variant_gene_counts_by_annotation"
ggsave(paste0(plot_title, '.png'), path = plots_dir, width=3.5, height=2)
ggsave(paste0(plot_title, '.pdf'), path = plots_dir, width=3.5, height=2)


###############################################################################
### Size changes of early-STOP genes ###
###############################################################################

# mydata$First_Annotation <- sub("&.*", "", mydata$Annotation)

nonsense <- mydata %>%
  # filter to stop gained variants
  subset(Annotation_simple == "Stop gained") %>%
  # create new columns from AA_pos column
  mutate(
    WT_AA_length = as.integer(sub(".*/", "", AA_pos)),
    nonsense_AA_length = as.integer(sub("/.*", "", AA_pos)),
    trunc_frac = 1 - (nonsense_AA_length / WT_AA_length)
  ) %>%
  # if a gene is in the list twice, filter to shorter nonsense product
  group_by(Gene_ID) %>%
  slice_min(nonsense_AA_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # set Gene_ID as a factor with levels sorted by truncation fraction
  mutate(Gene_ID = factor(Gene_ID, levels = Gene_ID[order(trunc_frac)])) %>%
  # make new column that is Gene_ID minus prefix (BV898_)
  mutate(Gene_ID_short = sub("BV898_", "", Gene_ID))
head(nonsense %>% select(Gene_ID, AA_pos, WT_AA_length, nonsense_AA_length, Uniprot_Protein_Annotation))

# ## ALL genes with Uniprot annotations
# ggplot(nonsense)+
#   geom_col(aes(x = Gene_ID, y = WT_AA_length), fill = "red", width = 0.8)+
#   geom_col(aes(x = Gene_ID, y = nonsense_AA_length), fill = "darkgrey", width = 0.8)+
#   geom_col(aes(x = Gene_ID, y = WT_AA_length), color = "black", fill = NA, width = 0.8)+
#   geom_text(
#     aes(
#       x = Gene_ID,
#       y = WT_AA_length * 1.03,   # a little past the bar end
#       label = Uniprot_Protein_Annotation
#     ),
#     hjust = 0,                    # left-align text relative to position
#     size = 3                      # adjust as needed
#   ) +
#   coord_flip()+
#   expand_limits(y = max(nonsense$WT_AA_length, na.rm = TRUE) * 2)+ # expand y axis for gene names
#   scale_y_continuous(
#     breaks = c(0,500,1000,1500,2000),
#     expand = c(0.01,0))+
#   labs(y = "Nonsense vs Wildtype Protein Length (AA)                                                ",
#        x = NULL,
#        title = "Uniprot Gene Annotations")
# plot_title <- "nonsense_gene_lengths_Uniprot_labels"
# ggsave(paste0(plot_title, '.png'), path = plots_dir, width=8, height=9)
# ggsave(paste0(plot_title, '.pdf'), path = plots_dir, width=8, height=9)

## Only genes with significant truncation, no annotations
nonsense_significant <- nonsense %>%
  filter(trunc_frac >= 0.25)  # keep only proteins truncated >= 25% (change if needed)
head(nonsense_significant %>% select(Gene_ID, AA_pos, WT_AA_length, nonsense_AA_length, trunc_frac))
ggplot(nonsense_significant)+
  geom_col(aes(x = Gene_ID_short, y = WT_AA_length), fill = "red", width = 0.8)+
  geom_col(aes(x = Gene_ID_short, y = nonsense_AA_length), fill = "darkgrey", width = 0.8)+
  geom_col(aes(x = Gene_ID_short, y = WT_AA_length), color = "black", fill = NA, width = 0.8)+
  coord_flip()+
  scale_y_continuous(
    breaks = c(0,500,1000,1500,2000))+
  labs(y = "Nonsense vs Wildtype Protein Length (AA)",
       x = "Locus ID (BV898_#)")
plot_title <- "nonsense_gene_lengths_significant_truncation_only"
ggsave(paste0(plot_title, '.png'), path = plots_dir, width=3.5, height=4.5)
ggsave(paste0(plot_title, '.pdf'), path = plots_dir, width=3.5, height=4.5)


###############################################################################
### Visualization of frameshift variants ###
###############################################################################

frameshift <- mydata %>%
  # filter to frameshift variants
  subset(Annotation_simple == "Frameshift") %>%
  # create new columns from AA_pos column
  mutate(
    WT_AA_length = as.integer(sub(".*/", "", AA_pos)),
    # nonsense length here means pre-frameshift length; just keeping same labels as above for ease
    nonsense_AA_length = as.integer(sub("/.*", "", AA_pos)),
    trunc_frac = 1 - (nonsense_AA_length / WT_AA_length)
  ) %>%
  # if a gene is in the list twice, filter to shorter nonsense product
  group_by(Gene_ID) %>%
  slice_min(nonsense_AA_length, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # set Gene_ID as a factor with levels sorted by WT_AA_length
  mutate(Gene_ID = factor(Gene_ID, levels = Gene_ID[order(trunc_frac)])) %>%
  # make new column that is Gene_ID minus prefix (BV898_)
  mutate(Gene_ID_short = sub("BV898_", "", Gene_ID))
head(frameshift %>% select(Gene_ID, Gene_ID_short, AA_pos, WT_AA_length, nonsense_AA_length, Uniprot_Protein_Annotation))

## Only genes with significantly early frameshift, no annotations
frameshift_significant <- frameshift %>%
  filter(trunc_frac >= 0.25)  # keep only proteins with frameshift >= 25% from C-terminus (change if needed)
head(frameshift_significant %>% select(Gene_ID, AA_pos, WT_AA_length, nonsense_AA_length, trunc_frac))
ggplot(frameshift_significant)+
  geom_col(aes(x = Gene_ID_short, y = WT_AA_length), fill = "blue", width = 0.8)+
  geom_col(aes(x = Gene_ID_short, y = nonsense_AA_length), fill = "darkgrey", width = 0.8)+
  geom_col(aes(x = Gene_ID_short, y = WT_AA_length), color = "black", fill = NA, width = 0.8)+
  coord_flip()+
  scale_y_continuous(
    breaks = c(0,500,1000,1500,2000))+
  labs(y = "Pre-frameshift vs Wildtype Protein Length (AA)",
       x = "Locus ID (BV898_#)")
plot_title <- "frameshift_gene_lengths_significant_proportion_only"
ggsave(paste0(plot_title, '.png'), path = plots_dir, width=3.5, height=4.5)
ggsave(paste0(plot_title, '.pdf'), path = plots_dir, width=3.5, height=4.5)


## save significant nonsense + frameshift genes to file
significant_truncations <- bind_rows(nonsense_significant, frameshift_significant) %>%
  arrange(Gene_ID) %>%
  select(Gene_ID, Annotation, Annotation_simple, AA_pos, WT_AA_length, nonsense_AA_length, trunc_frac, Uniprot_Protein_Annotation)
head(significant_truncations)
# save to file
write.csv(significant_truncations, "high_impact_variants_significant_truncations.csv", row.names = FALSE)

###############################################################################
### Analyze GO slim categories of select high-impact variants ###
###############################################################################

# pivot longer so each GO slim term is in its own row
goslim_long <- mydata %>%
  select(Gene_ID, GO_slim) %>%
  separate_rows(GO_slim, sep = ";")
head(goslim_long)
# count number of genes in each GO slim category
goslim_counts <- goslim_long %>%
  group_by(GO_slim) %>%
  summarise(n_genes = n_distinct(Gene_ID)) %>%
  arrange(desc(n_genes)) %>%
  # replace empty GO slim with "NA"
  mutate(GO_slim = ifelse(GO_slim == "" | is.na(GO_slim), NA, GO_slim))
head(goslim_counts)

# plot GO slim category counts
ggplot(goslim_counts)+
  geom_col(aes(x = reorder(GO_slim, n_genes), y = n_genes))+
  coord_flip()+
  labs(x = "GO slim category", y = "Number of genes",
       title = "Select high impact variants")
plot_title <- "select_high_impact_variants_GO_slim_category_counts"
ggsave(paste0(plot_title, '.png'), path = plots_dir, width=6, height=8)
ggsave(paste0(plot_title, '.pdf'), path = plots_dir, width=6, height=8)




