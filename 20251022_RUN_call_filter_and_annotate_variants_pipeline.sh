#!/bin/bash

set -euo pipefail
# Print helpful message if any command fails
trap 'echo "Error on line $LINENO: $BASH_COMMAND"; exit 1' ERR

module load bwa
module load samtools
module load bcftools
module load snpeff

###############################################
### User Settings ###
###############################################

script=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/20251010_wgs_call_filter_and_annotate_variants_pipeline.sh

# set per user
snpeff_dir=/nas/longleaf/home/adricoke/snpEff

# set per species: Hypsibius exemplaris
snpeff_database=Hypsibius_exemplaris_gca002082055v1
NCBI_ref=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/refs/Hypsibius_exemplaris/NCBI/GCA_002082055.1/GCA_002082055.1_nHd_3.1_genomic.fna
HiC_ref=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/refs/Hypsibius_exemplaris/HiC/nHd_3.1_HiC_5chromosomes.fa

###############################################
### Submit Jobs: Single-Tardigrade Sequencing vs NCBI Reference ###
###############################################

# cd /work/users/a/d/adricoke/single_tardigrade_seq_analysis/NCBI_genome_results

# for rep in 1 2 3 4
# do
#     if [ ! -d Hypsibius_exemplaris_rep${rep} ]; then
#         mkdir Hypsibius_exemplaris_rep${rep}
#     fi
#     cd Hypsibius_exemplaris_rep${rep}

#     # set per sample
#     reads=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/stwgs_rep${rep}.fastq.gz

#     sbatch --time=8:00:00 --mem=15g \
#     --job-name=NCBI_rep$rep --output=NCBI_rep${rep}.out \
#     --wrap="sh $script $NCBI_ref $reads $snpeff_database $snpeff_dir"

#     cd ..
# done

###############################################
### Submit Jobs: Single-Tardigrade vs Hi-C Genome ###
###############################################

# cd /work/users/a/d/adricoke/single_tardigrade_seq_analysis
# if [ ! -d new_HiC_genome_results ]; then
#     mkdir new_HiC_genome_results
# fi
# cd new_HiC_genome_results

# Use chromosome-level Hi-C genome as reference instead of NCBI assembly
# No SnpEff annotation available for this genome

# for rep in 1 2 3 4
# do
#     if [ ! -d Hypsibius_exemplaris_rep${rep} ]; then
#         mkdir Hypsibius_exemplaris_rep${rep}
#     fi
#     cd Hypsibius_exemplaris_rep${rep}

#     # set per sample
#     reads=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Arakawa2016_single_tardigrade_reads/stwgs_rep${rep}.fastq.gz

#     sbatch --time=8:00:00 --mem=15g \
#     --job-name=HiC_rep$rep --output=HiC_rep${rep}.out \
#     --wrap="sh $script $HiC_ref $reads"

#     cd ..
# done

###############################################
### Submit Jobs: Bulk-Tardigrade Sequencing vs NCBI Reference ###
###############################################

# cd /work/users/a/d/adricoke/single_tardigrade_seq_analysis/NCBI_genome_results

# ref=$NCBI_ref

# ### Boothby et al 2015 short reads ###

# READS_DIR=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Boothby2015_bulk_tardigrade_reads
# # imporant: assumes paired-end reads in format SAMPLE_1.fastq[.gz] and SAMPLE_2.fastq[.gz]
# declare -a SAMPLES=("SRR2986339" "SRR2986435" "SRR2986451")
# declare -a SAMPLE_NAMES=("bulk1" "bulk2" "bulk3")
# READ_EXT=".fastq"   # change to ".fastq.gz" if gzipped

# for idx in "${!SAMPLES[@]}"
# do
#     sample=${SAMPLES[$idx]}
#     sample_name=${SAMPLE_NAMES[$idx]}

#     if [ ! -d Hypsibius_exemplaris_${sample_name} ]; then
#         mkdir Hypsibius_exemplaris_${sample_name}
#     fi
#     cd Hypsibius_exemplaris_${sample_name}

#     reads_1="${READS_DIR}/${sample}_1${READ_EXT}"
#     reads_2="${READS_DIR}/${sample}_2${READ_EXT}"

#     sbatch --time=48:00:00 --mem=15g \
#     --job-name=${sample_name} --output=${sample_name}.out \
#     --wrap="sh $script $ref $reads_1 $reads_2 $snpeff_database $snpeff_dir"

#     cd ..
# done

###############################################
### Submit Jobs: Bulk-Tardigrade Sequencing vs HiC Reference ###
###############################################

cd /work/users/a/d/adricoke/single_tardigrade_seq_analysis/new_HiC_genome_results

ref=$HiC_ref

### Boothby et al 2015 short reads ###

READS_DIR=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Boothby2015_bulk_tardigrade_reads
# imporant: assumes paired-end reads in format SAMPLE_1.fastq[.gz] and SAMPLE_2.fastq[.gz]
declare -a SAMPLES=("SRR2986339" "SRR2986435" "SRR2986451")
declare -a SAMPLE_NAMES=("bulk1" "bulk2" "bulk3")
READ_EXT=".fastq"   # change to ".fastq.gz" if gzipped

for idx in "${!SAMPLES[@]}"
do
    sample=${SAMPLES[$idx]}
    sample_name=${SAMPLE_NAMES[$idx]}

    if [ ! -d Hypsibius_exemplaris_${sample_name} ]; then
        mkdir Hypsibius_exemplaris_${sample_name}
    fi
    cd Hypsibius_exemplaris_${sample_name}

    reads_1="${READS_DIR}/${sample}_1${READ_EXT}"
    reads_2="${READS_DIR}/${sample}_2${READ_EXT}"

    sbatch --time=48:00:00 --mem=15g \
    --job-name=bulk_HiC --output=${sample_name}.out \
    --wrap="sh $script $ref $reads_1 $reads_2"

    cd ..
done

# ### Yoshida et al 2017 Pacbio reads ###

# if [ ! -d Hypsibius_exemplaris_pacbio ]; then
#     mkdir Hypsibius_exemplaris_pacbio
# fi
# cd Hypsibius_exemplaris_pacbio

# # set per sample
# reads="/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Yoshida2017_bulk_tardigrade_PacBio_reads/SRR5179577_1.fastq"

# sbatch --time=48:00:00 --mem=15g \
# --job-name=pb --output=pb.out \
# --wrap="sh $script $ref $reads $snpeff_database $snpeff_dir"

# cd ..


###############################################
### Submit Jobs: Single-Animal Sequencing (NON-Tardigrade) ###
###############################################

# ### Drosophila melanogaster single-fly sample ###

# if [ ! -d Drosophila_melanogaster ]; then
#     mkdir Drosophila_melanogaster
# fi
# cd Drosophila_melanogaster

# # set per species
# ref=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/refs/Drosophila_melanogaster/NCBI/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
# snpeff_database=BDGP6.28.99

# # set per sample
# reads_1=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Adams2020_single_fly_reads/SRR10512945_1.fastq
# reads_2=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Adams2020_single_fly_reads/SRR10512945_2.fastq

# sbatch --time=48:00:00 --mem=15g \
# --job-name=fly --output=fly.out \
# --wrap="sh $script $ref $reads_1 $reads_2 $snpeff_database $snpeff_dir"

# cd ..

# ### Caenorhabditis elegans single-worm sample ###

# if [ ! -d Caenorhabditis_elegans ]; then
#     mkdir Caenorhabditis_elegans
# fi
# cd Caenorhabditis_elegans

# # set per species
# ref=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/refs/Caenorhabditis_elegans/NCBI/GCF_000002985.6/GCF_000002985.6_WBcel235_genomic.fna
# snpeff_database=WBcel235.99

# # set per sample
# reads_1=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Wang2024_single_worm_reads/cas2822_1.fq.gz
# reads_2=/work/users/a/d/adricoke/single_tardigrade_seq_analysis/data/Wang2024_single_worm_reads/cas2822_2.fq.gz

# sbatch --time=8:00:00 --mem=15g \
# --job-name=worm --output=worm.out \
# --wrap="sh $script $ref $reads_1 $reads_2 $snpeff_database $snpeff_dir"

# cd ..


###############################################
echo
echo "All jobs submitted."
echo