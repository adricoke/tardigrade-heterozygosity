#!/bin/bash
set -euo pipefail
trap 'echo "Error on line $LINENO: $BASH_COMMAND"; exit 1' ERR

module load bwa
module load samtools
module load gffread

###############################################
# USER SETTINGS
###############################################

### Use this script to map CDS sequences to the HiC genome assembly & NCBI reference assembly (for direct comparison)
    # saves cds_mask.csv files in each reference directory - list of covered positions in each CDS
    # assumes reference files are indexed with bwa index

# working directory
DIR=/work/users/a/d/adricoke/single_tardigrade_seq_analysis
cd $DIR

# NCBI reference genome and gene annotation file
NCBI_DIR=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/refs/Hypsibius_exemplaris/NCBI/GCA_002082055.1
NCBI_GFF=${NCBI_DIR}/genomic.gff
NCBI_REF=${NCBI_DIR}/GCA_002082055.1_nHd_3.1_genomic.fna

# HiC genome assembly
HIC_DIR=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/refs/Hypsibius_exemplaris/HiC
HIC_REF=${HIC_DIR}/nHd_3.1_HiC_5chromosomes.fa

###############################################
# GENERATE CDS FASTA FROM NCBI REFERENCE (ONCE)
###############################################
echo "[INFO] Generating cds.fa from NCBI GFF + reference..."
echo

if [[ ! -f "$NCBI_DIR/cds.fa" ]]; then
  echo "[INFO] Building cds.fa from GFF + reference..."
  echo
  gffread "$NCBI_GFF" -g "$NCBI_REF" -x "${NCBI_DIR}/cds.fa"
else
  echo "[INFO] cds.fa already exists; skipping."
  echo
fi

cds=${NCBI_DIR}/cds.fa

###############################################
# MAP CDS TO EACH REFERENCE GENOME
###############################################

# Array of reference directories and FASTA files
refs=(
  "$HIC_DIR:$HIC_REF"
  "$NCBI_DIR:$NCBI_REF"
)

for ref_entry in "${refs[@]}"; do
  # Split "dir:ref_fasta"
  ref_dir="${ref_entry%%:*}"
  ref="${ref_entry##*:}"

  echo "[INFO] Mapping CDS to reference genome:"
  echo "       Directory: $ref_dir"
  echo "       Reference: $ref"
  echo

  cd "$ref_dir"

  # map reads to ref with BWA algorithm, sort, and output as uncompressed bam file
  bwa mem "$ref" "$cds" | samtools sort -o cds_sorted.bam -

  # filter ambiguous mappings (MAPQ >= 20) & index filtered bam file
  samtools view -b -q 20 cds_sorted.bam > cds_filtered.bam
  samtools index cds_filtered.bam

  # Output a boolean CDS mask (chromosome, position, 1 = CDS-covered).
  # Only positions with CDS coverage are listed.
  samtools depth cds_filtered.bam \
    | awk '{print $1","$2",1"}' \
    > cds_mask.csv

done

###############################################
# DONE
###############################################
echo "[INFO] Pipeline complete."
echo "[INFO] Output cds_mask.csv files are located in each reference directory."
echo "[INFO] cds_mask.csv format: chromosome,position,CDS_coverage(1=yes). Only positions with CDS coverage are listed."
echo