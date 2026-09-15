#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=72:02:00
#SBATCH --output=%x-%j.out
#SBATCH --cpus-per-task=32

module load CBI

module load cellranger

cellranger count --id=run_count_LWHDRB2 \
    --fastqs=/c4/home/liw3/LWHDRB2/reads/GEX_cat \
    --sample=LWHDRB2_GEX \
    --transcriptome=/c4/home/liw3/ref_genome/refdata-gex-GRCm39-2024-A \
    --output-dir=/c4/home/liw3/LWHDRB2/processing_output/2_cellranger_output_cat \
    --create-bam=true \
    --localcores=32 \
    --localmem=128