# tardigrade-heterozygosity
Heterozygosity and variant analysis for the tardigrade *Hypsibius exemplaris*, as published in XXX.

# Code Descriptions

**Identify and analyze heterozygous variants:**

call_filter_and_annotate_variants_pipeline.sh
- Maps short genomic sequencing reads to reference with BWA algorithm
- Calls variants with bcftools
- Filters called variants by depth, quality, and whether they likely represent heterozygous sites
- Annotates variants with SnpEff
- Main outputs: heterozygous_variants.csv and all_variants.csv
- *Usage (per sample):*
   - Single-end reads: ./script.sh <reference.fa> <reads.fastq[.gz]> [snpeff_database snpeff_dir]
   - Paired-end reads: ./script.sh <reference.fa> <reads_1.fastq[.gz]> <reads_2.fastq[.gz]> [snpeff_database snpeff_dir]

map_cds_to_refs.sh
- Maps transcripts against a reference genome assembly for downstream CDS-only filtering
- Main output: cds_mask.csv

Compile_Heterozygosity_Data.R
- Compiles all_variants.csv files generated for multiple samples
- Computes the fraction of heterozygous sites out of total sites post-filtering, both in 100kb bins across putative chromosomes and total across the entire genome
- When applicable, applies CDS-only filter from cds_mask.csv
- Main outputs: heterozygosity_summary_data.csv and tardigrade_HiC_heterozygosity_all_samples_CDS_filter_FALSE_100kb_bins.csv

Plot_Heterozygosity_Data.R
- Reads in heterozygosity_summary_data.csv and tardigrade_HiC_heterozygosity_all_samples_CDS_filter_FALSE_100kb_bins.csv
- Makes plots of heterozygosity per sample and across putative chromosomes

Compile_Variant_Summary_Data.R
- Compiles heterozygous_variants.csv files generated for multiple samples
- Selects SnpEff annotation subfields of interest
- Generates summary file of variants counts per SnpEff annotation category: SnpEff_variant_counts.csv
- Filters variant list to those that are high impact and shared between all three bulk replicates and at least one single-animal replicate for downstream analyses: select_high_impact_variants.csv

Analyze_Tardigrade_High_Impact_Variants.R
- Reads in select_high_impact_variants.csv
- Makes plots to categorize and analyze select high-impact variants

**Analyze RNA expression of select variants:**

...

**Analyze results of PCR sequencing and inheritance tracking of select loci:**

...

# References

...