#!/bin/bash

#SBATCH -t 48:00:00
#SBATCH -c 6
#SBATCH --partition=biogipuzkoa-exclusive
#SBATCH --account=biogipuzkoa-exclusive
#SBATCH --mem=64G


module load Miniforge3
conda activate /scratch/mcarazo/envs/biotools

regions=("5s" "5.8s" "18s" "28s" "45s" "mt-12s" "mt-16s")

for bam in ../bams/*.bam; do
    #samtools index "$bam"
    echo "Processing $bam"
    sample=$(basename "$bam" .bam)
    for region in "${regions[@]}"; do
        bedtools genomecov -ibam "$bam" -dz > ${sample}.depth
        grep -w "$region" ${sample}.depth > ${sample}_${region}.depth
    done
    echo "Depth files created for $sample"
    echo "-----------------------------------"
done